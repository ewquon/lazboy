#!/usr/bin/env python
#
# GUI for semi-automatic SOWFA precursor set up
#
# Written by Eliot Quon (eliot.quon@nrel.gov), 2018-07-09
#
from __future__ import print_function
import os
import shutil
import yaml
import numpy as np

import tkinter as tk
from tkinter.scrolledtext import ScrolledText
try:
    # Python 2.7
    import tkFileDialog as filedialog
    import tkSimpleDialog as simpledialog
    import tkMessageBox as messagebox
except ImportError:
    # Python 3
    import tkinter.filedialog as filedialog
    import tkinter.simpledialog as simpledialog
    import tkinter.messagebox as messagebox

import template as tpl

DEBUG = True

pltstyle = 'seaborn-darkgrid'

# common files
soln_fields = ['U','T','p_rgh','k','kappat','nuSgs','Rwall']
constant_files = [
        'turbulenceProperties',
        'transportProperties',
        'LESProperties',
        'g',
        ] # write out custom forcingTable if needed
system_files = [
        'controlDict.1',
        'fvSolution',
        'fvSchemes',
        'decomposeParDict',
        'refineMeshDict.local',
        'changeDictionaryDict.updateBCs.cyclic'
        ] # write out custom setFieldABLDict
sampling_files = [
        'sliceDataInstantaneous',
        'boundaryDataPre'
        ]


#==============================================================================
#
# Definition of GUI and callback functions
#
class MainWindow(tk.Frame):
    """GUI for configuring the setUp file and generating a SOWFA
    precursor input deck
    """

    def __init__(self, master, default_template='neutral'):
        tk.Frame.__init__(self)

        self.master = master
        self.template = tk.StringVar(master, value=default_template)
        template_path = os.path.join(tpl.mypath,
                                     'simulation_templates',
                                     default_template+'.yaml')
        self.params = tpl.read_template(template_path)

        self.user_info()

        self.init_vars()
        self.init_window()
        self.restore_defaults()

        # resize window
        screen_w, screen_h = master.winfo_screenwidth(), root.winfo_screenheight()
        if DEBUG: print('Detected screen size:',screen_w,screen_h)
        master.update() # need to render window before gettin dimensions
        #master.attributes('-zoomed', True)
        widgets_w = self.scrollableframe.winfo_width() + self.scrollbar.winfo_width()
        master.geometry('{:d}x{:d}+0+0'.format(widgets_w,screen_h))



    def user_info(self):
        """Read or set up (first time) system information"""
        configpath = os.path.join(os.environ['HOME'],'.sowfa_defaults')
        try:
            with open(configpath,'r') as f:
                self.config = yaml.load(f)
        except IOError:
            self.config = dict()
            self.config['email'] = simpledialog.askstring(
                                        title='Default SOWFA configuration',
                                        prompt='Notification email',
                                        initialvalue='Richard.Astley@nrel.gov'
                                    )
            self.config['allocation'] = simpledialog.askstring(
                                            title='Default SOWFA configuration',
                                            prompt='Default allocation',
                                            initialvalue='windsim'
                                        )
            with open(configpath,'w') as f:
                yaml.dump(self.config, f, default_flow_style=False)
            print('Saved default SOWFA configuration info to '+configpath)
        if DEBUG: print(self.config)

            
    def init_vars(self):
        """Need to initialize some variables before we can create the
        widgets. However, most values associated with Entry boxes will
        be created on the fly.
        """
        self.lastrow = -1

        """decomposition controls"""
        self.decompTypeVar = tk.StringVar(self.master, value='simple')

        """general conditions"""
        self.coriolisVar = tk.IntVar(self.master, value=1)

        """atmosphere controls"""
        self.velocityInitTypeVar = tk.StringVar(self.master, value='geostrophic')
        self.temperatureInitTypeVar = tk.StringVar(self.master, value='simple')
        self.sourceTypeVar = tk.StringVar(self.master, value='constant')
        #self.p_rgh0 = 0.0  # Initial pressure (minus the hydrostatic variation and normalized by density) (m^2/s^2).
        #self.nuSgs0 = 0.0  # Initial SGS viscosity (m^2/s).
        #self.k0 = 0.1  # Initial SGS turbulent kinetic energy (m^2/s^2).
        #self.kappat0 = 0.0  # Initial SGS temperature diffusivity (m^2/s).
        self.idealProfileVar = tk.IntVar(self.master, value=0)
        self.momentumForcingVar = tk.IntVar(self.master, value=1)
        self.temperatureForcingVar = tk.IntVar(self.master, value=0)
 
        """surface controls"""
        self.surfaceBCTypeVar = tk.StringVar(self.master, value='fixed flux')
        #self.Rwall = [0.0,0.0,0.0,0.0,0.0,0.0]  # Initial wall shear stress (m^2/s^2).
        #self.kappa = 0.4  # von Karman constant.
 

    #--------------------------------------------------------------------------
    # TK wrappers
    #--------------------------------------------------------------------------

    def nextrow(self):
        self.lastrow += 1
        return self.lastrow

    def EntryRow(self,master,name,description,value=None,vcmd=None,extras=[],
                 **kwargs):
        if value is None:
            value = self.params[name]
        #setattr(self, name, tk.StringVar(value=str(value)))
        if vcmd is not None:
            kwargs['validate'] = 'focusout'
            #kwargs['validatecommand'] = (vcmd, name, '%W', '%P')
            kwargs['validatecommand'] = (vcmd,name)

        namestr = tk.Label(master, text=name)
        #entry = tk.Entry(master, textvariable=getattr(self,name), **kwargs)
        entry = tk.Entry(master, **kwargs)
        descstr = tk.Label(master, text=description)
        namestr.grid(row=self.nextrow(), column=0)
        entry.grid(row=self.lastrow, column=1)
        descstr.grid(row=self.lastrow, column=2, sticky=tk.W)
        for i,w in enumerate(extras):
            w.grid(row=self.lastrow, column=3+i)
        return entry

    def OptionMenuRow(self,master,name,description,args=[],variable=None,
                      **kwargs):
        namestr = tk.Label(master, text=name)
        entry = tk.OptionMenu(master, variable, *args, **kwargs)
        descstr = tk.Label(master, text=description)
        namestr.grid(row=self.nextrow(), column=0)
        entry.grid(row=self.lastrow, column=1, sticky='ew')
        descstr.grid(row=self.lastrow, column=2, sticky=tk.W)
        return entry

    def CheckboxRow(self,master,name,description,variable=None,**kwargs):
        namestr = tk.Label(master, text=name)
        entry = tk.Checkbutton(master, text='', variable=variable,
                               **kwargs)
        descstr = tk.Label(master, text=description)
        namestr.grid(row=self.nextrow(), column=0)
        entry.grid(row=self.lastrow, column=1, sticky='ew')
        descstr.grid(row=self.lastrow, column=2, sticky=tk.W)
        return entry

    #--------------------------------------------------------------------------
    # GUI setup
    #--------------------------------------------------------------------------

    def init_window(self):
        self.top = tk.Canvas(self.master)
        self.middle = tk.Canvas(self.master)
        self.bottom = tk.Canvas(self.master)
        self.top.pack(side='top', fill='x')
        self.middle.pack(side='top', fill='both', expand=True, pady=5)
        self.bottom.pack(side='top', fill='x', pady=5)

        # get available templates
        template_list, custom_config_list = tpl.get_templates()
        template_options = list(template_list.keys())
        if len(custom_config_list) > 0:
            template_options += ['---'] + list(custom_config_list.keys())
            for key,val in custom_config_list.items():
                template_list[key] = val
            self.template.set(key)
        self.template_list = template_list

        # top
        text = tk.Label(self.top, text='Canonical ABL template:')
        self.template_option = tk.OptionMenu(self.top,
                                             self.template,
                                             *template_options)
        self.restore_button = tk.Button(self.top,
                                        text='Restore defaults',
                                        command=self.restore_defaults)
        self.save_button = tk.Button(self.top,
                                     text='Save template',
                                     command=self.save_template)
        text.pack(side='left')
        self.template_option.pack(side='left', padx=20)
        self.restore_button.pack(side='left', padx=20)
        self.save_button.pack(side='left')

        # main part of window
        self.create_main_canvas()
        self.init_domain_controls()
        self.init_decomp_controls()
        self.init_general_controls()
        self.init_atmosphere_controls()
        self.init_surface_controls()
        #self.init_advanced_controls()

        # bottom
        self.generate_button = tk.Button(self.bottom,
                                         text='Generate case files!',
                                         fg='blue',
                                         command=self.generate)
        self.generate_button.pack(side='top', fill='x')


    def create_main_canvas(self):
        """Create a canvas to put all the controls in so we can
        implement a scrollbar.
        """
        self.canvas = tk.Canvas(self.middle)
        self.canvas.pack(side='left', fill='both', expand=True)
        self.scrollbar = tk.Scrollbar(self.middle, command=self.canvas.yview)
        self.scrollbar.pack(side='left', fill='y')
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        # update scrolling region after starting main loop
        def on_configure(event):
            self.canvas.configure(scrollregion=self.canvas.bbox('all'))
        self.canvas.bind('<Configure>', on_configure)

        # put frame in canvas
        # (this is the trick to make the scrollable frame work!)
        self.scrollableframe = tk.Frame(self.canvas)
        self.canvas.create_window((0,0), anchor='nw',
                                  window=self.scrollableframe)


    def init_domain_controls(self):
        #section = tk.LabelFrame(self.master, text='Domain')
        section = tk.LabelFrame(self.scrollableframe, text='Domain',
                                font='-weight bold')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)
        vcmd = self.register(self.update_grid_ext)
        self.xMin = self.EntryRow(section,
                                  'xMin', 'Minimum x-extent of domain [m]',
                                  vcmd=vcmd)
        self.yMin = self.EntryRow(section,
                                  'yMin', 'Minimum y-extent of domain [m]',
                                  vcmd=vcmd)
        self.zMin = self.EntryRow(section,
                                  'zMin', 'Minimum z-extent of domain [m]',
                                  vcmd=vcmd)
        self.xMax = self.EntryRow(section,
                                  'xMax', 'Maximum x-extent of domain [m]',
                                  vcmd=vcmd)
        self.yMax = self.EntryRow(section,
                                  'yMax', 'Maximum y-extent of domain [m]',
                                  vcmd=vcmd)
        self.zMax = self.EntryRow(section,
                                  'zMax', 'Maximum z-extent of domain [m]',
                                  vcmd=vcmd)

        self.dxText = tk.Label(section, text='dx=')
        self.dyText = tk.Label(section, text='dy=')
        self.dzText = tk.Label(section, text='dz=')
        self.NcellsText = tk.Label(section, text='Ncells=')

        vcmd2 = self.register(self.update_gridsize_decomp)
        self.nx = self.EntryRow(section, 'nx', 'Number of cells in x-direction',
                                vcmd=vcmd2,
                                extras=[self.dxText])
        self.ny = self.EntryRow(section, 'ny', 'Number of cells in y-direction',
                                vcmd=vcmd2,
                                extras=[self.dyText])
        self.nz = self.EntryRow(section, 'nz', 'Number of cells in z-direction',
                                vcmd=vcmd2,
                                extras=[self.dzText])
        self.NcellsText.grid(row=self.nextrow(), column=3)


    def init_decomp_controls(self):
        #section = tk.LabelFrame(self.master, text='Domain Decomposition')
        section = tk.LabelFrame(self.scrollableframe, text='Domain Decomposition',
                                font='-weight bold')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)
        vcmd = self.register(self.update_gridsize_decomp)
        self.PPN = self.EntryRow(section,
                                 'PPN', 'Processors per node',
                                 vcmd=vcmd)
        self.NnodesText = tk.Label(section, text='compute nodes')
        self.NnodesText.grid(row=self.lastrow, column=3)
        self.nCores = self.EntryRow(section,
                                    'nCores', 'Number of cores to use',
                                    vcmd=vcmd)
        self.cellsPerCoreText = tk.Label(section, text='cells/core')
        self.cellsPerCoreText.grid(row=self.lastrow, column=3)
        self.decompType = self.OptionMenuRow(section, 'decompType',
                                             'Decomposition algorithm',
                                             args=['simple','scotch'],
                                             variable=self.decompTypeVar,
                                             command=self.update_decomp)
        self.decompOrder = self.EntryRow(section,
                                         'decompOrder',
                                         'Order of decomposition in x,y,z-directions',
                                         )


    def init_general_controls(self):
        #section = tk.LabelFrame(self.master, text='General Settings')
        section = tk.LabelFrame(self.scrollableframe, text='General Settings',
                                font='-weight bold')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)
        self.TRef = self.EntryRow(section,
                                  'TRef','Reference potential temperature (K)')
        self.coriolis = self.CheckboxRow(section,'coriolis','Coriolis forces',
                                         variable=self.coriolisVar,
                                         command=self.update_coriolis)
        self.latitude = self.EntryRow(section,
                                      'latitude',
                                      'Latitude on the Earth of the site (deg)')
        self.EarthPeriod = self.EntryRow(section,
                                         'EarthPeriod',
                                         'Earth\'s rotation period (hr)')


    def init_atmosphere_controls(self):
        #section = tk.LabelFrame(self.master, text='Atmospheric Conditions')
        section = tk.LabelFrame(self.scrollableframe, text='Atmospheric Conditions',
                                font='-weight bold')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)

        header = tk.Label(section,
                          text='Initial conditions', font='-slant italic')
        header.grid(row=self.nextrow(), column=0, sticky=tk.W)
        self.plot_init_button = tk.Button(section,
                                          text='PLOT',
                                          command=self.plot_initial_conditions)
        self.plot_init_button.grid(row=self.lastrow, column=3, sticky=tk.E)

        self.velocityInitType = self.OptionMenuRow(section, 'velocityInitType',
                                                   'How to initialize the base velocity profile',
                                                   args=['geostrophic','log','table'],
                                                   variable=self.velocityInitTypeVar,
                                                   command=self.update_velocity_init)
        vcmd = self.register(self.update_profiles)
        self.U0Mag = self.EntryRow(section,
                                   'U0Mag',
                                   'Initial condition for wind speed (m/s).',
                                   vcmd=vcmd)
        self.dir = self.EntryRow(section,
                                 'dir',
                                 'Initial condition for wind direction (deg).',
                                 vcmd=vcmd)
        self.windHeight = self.EntryRow(section,
                                        'windHeight',
                                        'Height at which to drive mean wind to U0Mag/dir (m).',
                                        vcmd=vcmd)
        #self.p_rgh0 = 0.0  # Initial pressure (minus the hydrostatic variation and normalized by density) (m^2/s^2).
        #self.nuSgs0 = 0.0  # Initial SGS viscosity (m^2/s).
        #self.k0 = 0.1  # Initial SGS turbulent kinetic energy (m^2/s^2).
        #self.kappat0 = 0.0  # Initial SGS temperature diffusivity (m^2/s).

        self.temperatureInitType = self.OptionMenuRow(section, 'temperatureInitType',
                                                   'How to initialize the base temperature profile',
                                                   args=['simple','table'],
                                                   variable=self.temperatureInitTypeVar,
                                                   command=self.update_temperature_init)
        self.TGradUpper = self.EntryRow(section,
                                        'TGradUpper',
                                        'Potential temperature gradient above the strong inversion (K/m).',
                                        vcmd=vcmd)
        self.zInversion = self.EntryRow(section,
                                        'zInversion',
                                        'Height of the middle of the initial strong capping inversion (m).',
                                        vcmd=vcmd)
        self.inversionWidth = self.EntryRow(section,
                                            'inversionWidth',
                                            'Vertical width of the intial strong capping inversion (m).',
                                            vcmd=vcmd)
        self.TBottom = self.EntryRow(section,
                                     'TBottom',
                                     'Initial potential temperature at bottom of strong capping inversion (K).',
                                     vcmd=vcmd)
        self.TTop = self.EntryRow(section,
                                  'TTop',
                                  'Initial potential temperature at top of strong capping inversion (K).',
                                  vcmd=vcmd)
        self.TgradText = tk.Label(section, text='Tgrad = ')
        self.TgradText.grid(row=self.nextrow(), column=2)

        header = tk.Label(section,
                          text='Background conditions', font='-slant italic')
        header.grid(row=self.nextrow(), column=0, sticky=tk.W)
        self.plot_bkgnd_button = tk.Button(section,
                                           text='PLOT',
                                           command=self.plot_background_conditions)
        #self.plot_bkgnd_button.grid(row=self.lastrow, column=3, sticky=tk.E)

        self.sourceType = self.OptionMenuRow(section, 'sourceType',
                                             'Computed source type',
                                             args=['constant','profile'],
                                             variable=self.sourceTypeVar,
                                             command=self.update_source_type)
        self.idealProfile = self.CheckboxRow(section,'idealProfile',
                                             'Compute velocity profile from U0Mag, dir, windHeight, alpha, and veer',
                                             variable=self.idealProfileVar,
                                             command=self.update_ideal_profile)
        self.momentumForcing = self.CheckboxRow(section,'momentumForcing',
                                                'Adjust source terms to match U0Mag at windHeight or specified velocity profile',
                                                variable=self.momentumForcingVar,
                                                state='disabled')
        self.temperatureForcing = self.CheckboxRow(section,'temperatureForcing',
                                                   'Adjust source terms to match specified temperature profile',
                                                   variable=self.temperatureForcingVar)
        vcmd = self.register(self.update_profiles)
        self.alpha = self.EntryRow(section, 'alpha', 'Shear exponent',
                                   vcmd=vcmd)
        self.veer = self.EntryRow(section, 'veer', 'Veer over windHeight (deg)',
                                  vcmd=vcmd)

        """Optional profile table for finer control over the background
        and initial conditions. Note, I'm assuming that if a profile is
        specified for the background conditions, the same profile should
        be used to initialize the solution with setFieldsABL. Also, the
        background conditions are assumed to be stationary.
        """
        profileTableText = 'profileTable\n' + \
                '(z, U, V, theta)\n\n' + \
                '***uncheck idealProfile\nif manually edited***'
        namestr = tk.Label(section, text=profileTableText)
        namestr.grid(row=self.nextrow(), column=0)
        text = ScrolledText(section, borderwidth=1)
        text.grid(row=self.lastrow, column=1, columnspan=2, sticky='ew')
        self.profileTable = text
        self.plot_bkgnd_button.grid(row=self.lastrow, column=3, sticky=tk.E)


    def init_surface_controls(self):
        #section = tk.LabelFrame(self.master, text='Surface Conditions')
        section = tk.LabelFrame(self.scrollableframe, text='Surface Conditions',
                                font='-weight bold')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)
        self.surfaceBCType = self.OptionMenuRow(
                section,
                'surfaceBCType',
                'Surface conditions (flux for unstable, heating rate for stable)',
                args=['fixed flux','fixed heating rate'],
                variable=self.surfaceBCTypeVar,
                command=self.update_surface_bc)
        vcmd = self.register(self.update_surface_bc)
        self.z0 = self.EntryRow(section, 'z0',
                                'Surface roughness, for loglaw init and Rwall (m)',
                                vcmd=vcmd)
        self.fluxInfoText = tk.Label(section, text='= W/m^2')
        self.qwall = self.EntryRow(section, 'qwall',
                                   'Temperature flux, negative is into domain (K-m/s)',
                                   vcmd=vcmd,
                                   extras=[self.fluxInfoText])
        self.heatingInfoText = tk.Label(section, text='heating')
        self.heatingRate = self.EntryRow(section, 'heatingRate',
                                         'Surface temperature change rate, negative is cooling (K/s)',
                                         vcmd=vcmd,
                                         extras=[self.heatingInfoText])

#    def init_advanced_controls(self):
#        #section = tk.LabelFrame(self.master, text='Advanced Controls')
#        section = tk.LabelFrame(self.scrollableframe, text='Advanced Controls',
#                                font='-weight bold')
#        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)

    #--------------------------------------------------------------------------
    # validation routines (callback functions)
    #--------------------------------------------------------------------------

    def update_grid_ext(self,name):
        newval = getattr(self,name).get()
        try:
            fval = float(newval)
        except ValueError:
            self.bell()
            #getattr(self,name).set(self.params[name]) # THIS BREAKS THE VALIDATION
            getattr(self,name).delete(0,tk.END)
            getattr(self,name).insert(0,self.params[name])
            return False # this is meaningless for 'focusout' validation
        else:
            self.params[name] = fval
            self.calc_grid_res()
            return True

    def update_gridsize_decomp(self,name):
        newval = getattr(self,name).get()
        try:
            ival = int(newval)
            assert(ival > 0)
        except (ValueError, AssertionError):
            self.bell()
            getattr(self,name).delete(0,tk.END)
            getattr(self,name).insert(0,self.params[name])
            return False
        else:
            if name=='nCores':
                ppn = int(self.PPN.get())
                if not (ival % ppn == 0):
                    ival = int(np.round(float(ival) / ppn) * ppn)
                    getattr(self,name).delete(0,tk.END)
                    getattr(self,name).insert(0,ival)
            self.params[name] = ival
            self.calc_grid_res()
            return True

    def calc_grid_res(self):
        try:
            x0 = float(self.xMin.get())
            y0 = float(self.yMin.get())
            z0 = float(self.zMin.get())
            x1 = float(self.xMax.get())
            y1 = float(self.yMax.get())
            z1 = float(self.zMax.get())
            nx = int(self.nx.get())
            ny = int(self.ny.get())
            nz = int(self.nz.get())
        except ValueError:
            pass
        else:
            self.dx = (x1-x0) / nx
            self.dy = (y1-y0) / ny
            self.dz = (z1-z0) / nz
            self.zsurf = z0
            self.dxText['text'] = 'dx = {} m'.format(self.dx)
            self.dyText['text'] = 'dy = {} m'.format(self.dy)
            self.dzText['text'] = 'dz = {} m'.format(self.dz)
            self.NcellsText['text'] = 'Ncells = {:g}'.format(nx*ny*nz)

            self.coresNeeded = int(self.nCores.get())/int(self.PPN.get())
            self.NnodesText['text'] = '{} compute nodes'.format(self.coresNeeded)

            self.avgCellsPerCore = int(nx*ny*nz/float(self.nCores.get()))
            self.cellsPerCoreText['text'] = 'average {} cells/core'.format(self.avgCellsPerCore)


    def update_profiles(self,name=None):
        """Update basic atmospheric properties, including the ideal
        temperature profile for ICs. If an ideal profile is specified,
        then horizontal wind components described by shear and veer
        parameters is also calculated.
        """
        zref = float(self.windHeight.get())
        zinv = float(self.zInversion.get())
        Tbot = float(self.TBottom.get())
        Ttop = float(self.TTop.get())
        width = float(self.inversionWidth.get())
        Tgrad = float(self.TGradUpper.get())
        wdir = float(self.dir.get())

        self.Tgrad_strong = 100 * (Ttop-Tbot) / width
        self.TgradText['text'] = 'strong inversion gradient: {:.1f} K/(100 m)'.format(self.Tgrad_strong)

        # for plots
        nz = int(self.nz.get())
        self.z = self.zsurf + self.dz/2 + self.dz*np.arange(nz)

        # for user-specified init/background profile
        listlist = self._text_to_listlist(self.profileTable.get('1.0',tk.END))
        profiledata = np.array(listlist)

        # temperature profile
        temperatureInit = self.temperatureInitTypeVar.get()
        if (temperatureInit == 'simple') or (len(profiledata) == 0):
            self.z_T = self.z
            zbot = zinv - width/2
            ztop = zinv + width/2
            strong_layer = (self.z >= zbot) & (self.z <= ztop)
            weak_layer = (self.z > ztop)
            self.T = Tbot * np.ones(self.z.shape)
            self.T[strong_layer] = Tbot + (self.z[strong_layer] - zbot) * (Ttop-Tbot)/width
            self.T[weak_layer] = Ttop + (self.z[weak_layer] - ztop) * Tgrad
        elif temperatureInit == 'table':
            self.z_T = profiledata[:,0]
            self.T = profiledata[:,3]
        else:
            raise ValueError('Unknown temperature initialization: {:s}'.format(temperatureInit))

        if self.idealProfileVar.get():
            # calculate ideal wind profile
            self.z_U = self.z
            Uref = float(self.U0Mag.get())
            alpha = float(self.alpha.get())
            veer = float(self.veer.get())

            self.WS = Uref * (self.z / zref)**alpha
            self.WD = wdir + veer/zref * (self.z - zref)
            freeatm = self.z >= zinv
            self.WD[freeatm] = self.WD[freeatm][0]

            ang = 1.5*np.pi - np.pi/180*self.WD
            self.U = self.WS * np.cos(ang)
            self.V = self.WS * np.sin(ang)

            self._update_profile_text()

        else:
            # get input wind profile
            self.z_U = profiledata[:,0]
            self.U = profiledata[:,1]
            self.V = profiledata[:,2]
            self.WS = np.sqrt(self.U**2 + self.V**2)
            self.WD = 180./np.pi * np.arctan2(-self.U, -self.V)
            self.WD[self.WD < 0] += 360.

        # update wind direction
        sourceType = self.sourceTypeVar.get()
        if sourceType=='profile':
            # if given profile, interpolate wind direction at hub height
            wdir = np.interp(zref, self.z_U, self.WD)
        ang = 1.5*np.pi - wdir*np.pi/180.
        self.wdir_vec = np.array([np.cos(ang), np.sin(ang), 0])
        #if (wdir >= 337.5) or (wdir < 22.5):
        #    self.wdirname = 'north'
        #elif (wdir >= 22.5) and (wdir < 67.5):
        #    self.wdirname = 'northeast'
        #elif (wdir >= 67.5) and (wdir < 112.5):
        #    self.wdirname = 'east'
        #elif (wdir >= 112.5) and (wdir < 157.5):
        #    self.wdirname = 'southeast'
        #elif (wdir >= 157.5) and (wdir < 202.5):
        #    self.wdirname = 'south'
        #elif (wdir >= 202.5) and (wdir < 247.5):
        #    self.wdirname = 'southwest'
        #elif (wdir >= 247.5) and (wdir < 292.5):
        #    self.wdirname = 'west'
        #elif (wdir >= 292.5) and (wdir < 337.5):
        #    self.wdirname = 'northwest'
        # These are not true compass sectors, but rather, directions
        # corresponding to SOWFA inflow
        if wdir == 0.0:
            self.wdirname = 'north'
        elif (wdir > 0.0) and (wdir < 90.0):
            self.wdirname = 'northeast'
        elif wdir == 90.0:
            self.wdirname = 'east'
        elif (wdir > 90.0) and (wdir < 180.0):
            self.wdirname = 'southeast'
        elif wdir == 180.0:
            self.wdirname = 'south'
        elif (wdir > 180.0) and (wdir < 270.0):
            self.wdirname = 'southwest'
        elif wdir == 270.0:
            self.wdirname = 'west'
        elif (wdir > 270.0) and (wdir < 360.0):
            self.wdirname = 'northwest'
        else:
            raise ValueError('Unexpected wind direction: {:f}'.format(wdir))
        if DEBUG:
            print('Wind direction (z={:g} m) : {:g}, SOWFA "{:s}erly" wind'.format(zref, wdir,self.wdirname))

        return True # so the widget isn't disabled


    def _update_profile_text(self):
        self.profileTable.delete('1.0', tk.END)
        #self.profileTable.insert(tk.END, '//   z       U       V       T\n')
        for zi,Ui,Vi,Ti in zip(self.z,self.U,self.V,self.T):
            rowdata = '{}  {}  {}  {}\n'.format(zi,Ui,Vi,Ti)
            self.profileTable.insert(tk.END, rowdata)


    def _disable_profile_table(self):
        if not ( (self.sourceTypeVar.get() == 'profile') or
                 (self.velocityInitTypeVar.get() == 'table') or
                 (self.temperatureInitTypeVar.get() == 'table') ):
            self.profileTable.config(state='disabled',fg='light grey')
            self.plot_bkgnd_button.config(state='disabled')


    def update_decomp(self,decompType):
        if decompType=='simple':
            self.decompOrder.config(state='normal')
        elif decompType=='scotch':
            self.decompOrder.config(state='disabled')
        else:
            raise ValueError('Unexpected decompType: '+decompType)


    def update_coriolis(self):
        coriolis = self.coriolisVar.get()
        if coriolis==1:
            self.EarthPeriod.config(state='normal')
            self.latitude.config(state='normal')
        else:
            self.EarthPeriod.config(state='disabled')
            self.latitude.config(state='disabled')


    def update_ideal_profile(self):
        sourceType = self.sourceTypeVar.get()
        ideal = self.idealProfileVar.get()
        if (sourceType=='profile') and (ideal==1):
            self.alpha.config(state='normal')
            self.veer.config(state='normal')
        else:
            self.alpha.config(state='disabled')
            self.veer.config(state='disabled')

    def update_source_type(self,sourceType):
        if sourceType=='constant':
            self._disable_profile_table()
            #self.momentumForcing.config(state='disabled')
            self.temperatureForcing.config(state='disabled')
            self.idealProfile.config(state='disabled')
            self.alpha.config(state='disabled')
            self.veer.config(state='disabled')
        elif sourceType=='profile':
            self.profileTable.config(state='normal',fg='black')
            self.plot_bkgnd_button.config(state='normal')
            #self.momentumForcing.config(state='normal')
            self.temperatureForcing.config(state='normal')
            self.idealProfile.config(state='normal')
            self.alpha.config(state='normal')
            self.veer.config(state='normal')
        else:
            raise ValueError('Unexpected source type: '+sourceType)
        self.update_profiles()


    def update_surface_bc(self,name=None):
        R = 287. # [J/(kg-K)]
        pref = 101300. # standard atmospheric pressure at sea level [Pa]
        rhoref = pref / R / float(self.TRef.get())
        Cp = 3.5*R # gamma/(gamma-1) * R
        qwall = float(self.qwall.get())
        self.heating = float(self.heatingRate.get())

        self.flux = qwall * rhoref * Cp
        fluxdesc = ''
        if self.flux < 0:
            fluxdesc = 'HEATING'
        elif self.flux > 0:
            fluxdesc = 'COOLING'
        self.fluxInfoText['text'] = '= {} W/m^2 {}'.format(np.abs(self.flux),
                                                           fluxdesc)

        heatingdesc = ''
        if self.heating > 0:
            heatingdesc = 'HEATING'
        elif self.heating < 0:
            heatingdesc = 'COOLING'
        self.heatingInfoText['text'] = '= {:.4f} K/hr {}'.format(3600*np.abs(self.heating),
                                                             heatingdesc)

        surfaceBCType = self.surfaceBCTypeVar.get()
        if surfaceBCType == 'fixed flux':
            self.qwall.config(state='normal')
            self.heatingRate.config(state='disabled')
        elif surfaceBCType == 'fixed heating rate':
            self.qwall.config(state='disabled')
            self.heatingRate.config(state='normal')
        else:
            raise ValueError('Unexpected surfaceBCType: '+surfaceBCType)

        return True # so the widget isn't disabled


    def update_velocity_init(self,initType):
        if initType == 'table':
            self.plot_bkgnd_button.config(state='normal')
            self.profileTable.config(state='normal',fg='black')
            self.update_profiles()
        else:
            self._disable_profile_table()

    def update_temperature_init(self,initType):
        if initType == 'simple':
            #self.TGradUpper.config(state='normal')
            self.zInversion.config(state='normal')
            self.inversionWidth.config(state='normal')
            self.TBottom.config(state='normal')
            self.TTop.config(state='normal')
            self._disable_profile_table()
        elif initType == 'table':
            #self.TGradUpper.config(state='disabled')
            self.zInversion.config(state='disabled')
            self.inversionWidth.config(state='disabled')
            self.TBottom.config(state='disabled')
            self.TTop.config(state='disabled')
            self.plot_bkgnd_button.config(state='normal')
            self.profileTable.config(state='normal',fg='black')
            self.update_profiles()
        else:
            raise ValueError('Unexpected temperatureInitType: '+initType)


    #--------------------------------------------------------------------------
    # actions
    #--------------------------------------------------------------------------

    def _text_to_list(self,line,dtype=float):
        """Parse "[1 2 3]" into a list."""
        L = [ dtype(val) for val in line.strip().strip('[]').split(',') ]
        return L

    def _text_to_listlist(self,text):
        """Parse
        0.0 0.0 0.0 273.121917725
        5.0 1.82671529027 1.98969371553 282.524156211
        10.0 2.52983624172 2.77228242751 283.156504284
        15.0 3.01282071998 3.35521702288 283.64108812
        ...
        into a list of lists, for conversion to 2-D np.ndarray.
        """
        lines = text.split('\n')
        #listlist = [ [ float(val) for val in line.strip().strip('()').split() ]
        listlist = [ [ float(val) for val in line.strip().split() ]
                    for line in lines if not line.strip()=='' ]
        return listlist

    def _listlist_to_profiles(self,listlist):
        if (listlist is not None) and (len(listlist) > 1):
            data = np.array(listlist)
            self.z = data[:,0]
            self.U = data[:,1]
            self.V = data[:,2]
            self.T = data[:,3]
            self._update_profile_text()


    def restore_defaults(self):
        """It's kind of a PITA, but binding tk.StringVar to the widgets
        and trying to udpate it in the validate command breaks the
        callback routine. So instead of being able to just call 
            tk.StringVar.set(value),
        we instead have to call
            tk.Entry.delete(0, END)
            tk.Entry.insert(0, value).
        """
        name = self.template.get()
        print('Restoring defaults for '+name)
        params = tpl.read_template(self.template_list[name])
        if DEBUG: print(params)
        self.params = params

        for name, val in params.items():
            try:
                widget = getattr(self, name)
            except AttributeError:
                print('Note: widget "'+name+'" does not exist')
            else:
                widget_class = widget.winfo_class()
                if widget_class == 'Entry':
                    if DEBUG: print('Setting entry "'+name+'"')
                    widget.delete(0, tk.END)
                    widget.insert(0, str(val))
                elif widget_class == 'Text':
                    if DEBUG:
                        print('Updating profiles and setting "'+name+'"')
                        print('  list of lists: '+str(val))
                    self._listlist_to_profiles(val)
                else:
                    # we have a control variable instead of widget
                    if DEBUG: print('Setting variable for "{}" ({})'.format(name,widget_class))
                    wvar = getattr(self, name+'Var')
                    wvar.set(val)
        
        """update active text widgets"""
        # TODO: update these functions with defaults so we don't need to call *.get()
        self.calc_grid_res()
        #self.update_profiles()
        self.update_source_type(self.sourceTypeVar.get()) # calls update_profiles()
        self.update_decomp(self.decompTypeVar.get())
        self.update_coriolis()
        self.update_velocity_init(self.velocityInitTypeVar.get())
        self.update_temperature_init(self.temperatureInitTypeVar.get())
        self.update_ideal_profile()
        self.update_surface_bc()


    def get_all_params(self):
        """Called before saving template and generating case files to
        update the 'params' dictionary.
        """
        for name, val in self.params.items():
            dtype = type(val)
            try:
                widget = getattr(self, name)
            except AttributeError:
                pass
            else:
                if widget.winfo_class() == 'Entry':
                    if DEBUG: print('Updated variable "{}" {} from Entry'.format(name,str(dtype)))
                    wvar = widget.get()
                    if dtype is list:
                        dtype = type(val[0])
                        vals = self._text_to_list(wvar,dtype=dtype)
                        if DEBUG: print('  reconstructed list: {}'.format(vals))
                        self.params[name] = vals
                    else:
                        self.params[name] = dtype(wvar)
                elif widget.winfo_class() == 'Text':
                    if DEBUG: print('Updated variable "{}" {} from Text'.format(name,str(dtype)))
                    text = widget.get('1.0',tk.END)
                    #listlist = self._text_to_listlist(text)
                    #if DEBUG: print(' reconstructed list of lists: {}'.format(listlist))
                    #self.params[name] = listlist
                    self.params[name] = text.strip()
                else:
                    # need to get from control variable instead of widget
                    if DEBUG:
                        print('Updated variable "{}" {} from {}'.format(
                            name, str(dtype), widget.winfo_class()) )
                    wvar = getattr(self, name+'Var').get()
                    self.params[name] = dtype(wvar)


    def save_template(self,fpath=None):
        """Save yaml file which should completely describe the
        simulation
        """
        if fpath is None:
            self.get_all_params()
            fpath = filedialog.asksaveasfilename(
                        initialdir=os.path.join(tpl.mypath, 'simulation_templates'),
                        title='Save precursor configuration',
                        filetypes=(('yaml files','*.yaml'),('all files','*,*')),
                    )
        params_copy = self.params.copy()

        # clean up variables for output as needed
        listlist = self._text_to_listlist(params_copy['profileTable'])
        if DEBUG: print(' reconstructed list of lists: {}'.format(listlist))
        params_copy['profileTable'] = listlist

        # now write it out
        if not fpath == '':
            with open(fpath,'w') as f:
                f.write(tpl.yaml_template.format(**params_copy))
        print('Wrote out '+fpath)


    def generate(self):
        """Create the OpenFOAM simulation structure. Called when user
        clicks the 'Generate case files!' button.
        """
        self.get_all_params()
        for key,val in self.config.items():
            self.params[key] = val
        self.params['nNodes'] = int(self.coresNeeded)

        self.check_sanity()

        dpath = filedialog.askdirectory(
                    initialdir=os.getcwd(),
                    title='Specify new simulation directory',
                )

        """Create run directory structure here"""
        if not dpath == '':
            if not os.path.isdir(dpath):
                os.makedirs(dpath)
            casename = os.path.split(dpath)[-1]
            self.params['casename'] = simpledialog.askstring(
                                            title='runscript',
                                            prompt='PBS job name',
                                            initialvalue=casename
                                        )
            srcdir = os.path.join(tpl.mypath, 'simulation_templates')
            print('Generating case files in '+dpath)
            print('  from files in '+srcdir)
            self.save_template(fpath=os.path.join(dpath,'setup.yaml'))

            # fix for booleans (OpenFOAM expects lowercase)
            params_copy = self.params.copy()
            for key,val in params_copy.items():
                if isinstance(val, bool):
                    params_copy[key] = str(val).lower()

            # additional calculations for substitution into sowfa templates
            off = 0.01
            params_copy['halfx'] = (self.params['xMin']+self.params['xMax']) / 2 + off
            params_copy['halfy'] = (self.params['yMin']+self.params['yMax']) / 2 + off
            params_copy['halfz'] = (self.params['zMin']+self.params['zMax']) / 2 + off
            params_copy['zhub'] = self.params['windHeight'] + off

            norm = np.cross([0,0,1], self.wdir_vec)
            norm_components = [ round(val,8) for val in norm ]
            params_copy['vertical_plane_normal'] = '{:g} {:g} {:g}'.format(*norm_components)

            # select boundary patches
            inlet_patches = 'surfaces\n          (\n'
            if len(self.wdirname) <= 5:
                # single inlet
                patches = [self.wdirname]
            else:
                # inflow not aligned with a single patch
                patches = [self.wdirname[:5], self.wdirname[5:]]
            print('Inflow patches: {}'.format(patches))
            for patch in patches:
                inlet_patches += """              {name}
              {{
                  type         patch;
                  patches      ({name});
                  triangulate  false;
              }}\n""".format(name=patch)
            inlet_patches += '          );'
            params_copy['inlet_patches'] = inlet_patches

            # convert profile table to OpenFOAM format
            data = ''
            for line in params_copy['profileTable'].split('\n'):
                if not line == '':
                    data += '(' + line + ')\n'
            params_copy['profileTable'] = data.strip()

            # write setUp file
            tpl.copy_and_update(srcdir,dpath,
                                'setUp',
                                params_copy)

            # initial conditions
            os.makedirs(os.path.join(dpath,'0.original'))
            for fname in soln_fields:
                shutil.copy2(os.path.join(srcdir,'0.original',fname),
                             os.path.join(dpath,'0.original'))
            surf_bc = self.surfaceBCTypeVar.get().replace(' ','')
            shutil.copy2(os.path.join(srcdir,'0.original','qwall.'+surf_bc),
                         os.path.join(dpath,'0.original','qwall'))

            # constant files
            os.makedirs(os.path.join(dpath,'constant','polyMesh'))
            for fname in constant_files:
                shutil.copy2(os.path.join(srcdir,'constant',fname),
                             os.path.join(dpath,'constant'))
            shutil.copy2(os.path.join(srcdir,'constant','polyMesh','blockMeshDict'),
                         os.path.join(dpath,'constant','polyMesh'))
            source_type = self.sourceTypeVar.get().replace(' ','')
            if source_type == 'constant':
                shutil.copy2(os.path.join(srcdir,'constant','ABLProperties.'+source_type),
                             os.path.join(dpath,'constant','ABLProperties'))
            elif source_type == 'profile':
                # setup sources
                if self.momentumForcingVar.get():
                    params_copy['momentumSourceType'] = 'computed'
                    params_copy['momentumForcing'] = '#include "momentumForcingTable"'
                    # write out momentum forcing table
                    fpath = os.path.join(dpath,'constant','momentumForcingTable')
                    xmom_sources = ' '.join([str(u) for u in self.U])
                    ymom_sources = ' '.join([str(v) for v in self.V])
                    zmom_sources = ' '.join(len(self.z_U)*['0.0'])
                    assert(len(self.U) == len(self.V) == len(self.z_U))
                    with open(fpath,'w') as f:
                        f.write('sourceHeightsMomentum\n(\n\t')
                        f.write('\n\t'.join([str(z) for z in self.z_U]))
                        f.write('\n);\n\n')
                        f.write('sourceTableMomentumX\n(\n')
                        f.write('\t(0.0 '+xmom_sources+')\n')
                        f.write('\t(90000.0 '+xmom_sources+')\n')
                        f.write(');\n\n')
                        f.write('sourceTableMomentumY\n(\n')
                        f.write('\t(0.0 '+ymom_sources+')\n')
                        f.write('\t(90000.0 '+ymom_sources+')\n')
                        f.write(');\n\n')
                        f.write('sourceTableMomentumZ\n(\n')
                        f.write('\t(0.0 '+zmom_sources+')\n')
                        f.write('\t(90000.0 '+zmom_sources+')\n')
                        f.write(');\n\n')
                    print('Wrote out '+fpath)
                else:
                    params_copy['momentumSourceType'] = 'given'
                    self._alert('source_type=="profile", momentumSourceType=="given" not handled')

                if self.temperatureForcingVar.get():
                    params_copy['temperatureSourceType'] = 'computed'
                    params_copy['temperatureForcing'] = '#include "temperatureForcingTable"'
                    # write out temperature forcing table
                    fpath = os.path.join(dpath,'constant','temperatureForcingTable')
                    T_sources = ' '.join([str(T) for T in self.T])
                    assert(len(self.T) == len(self.z_T))
                    with open(fpath,'w') as f:
                        f.write('sourceHeightsTemperature\n(\n\t')
                        f.write('\n\t'.join([str(z) for z in self.z_T]))
                        f.write('\n);\n\n')
                        f.write('sourceTableTemperature\n(\n')
                        f.write('\t(0.0 '+T_sources+')\n')
                        f.write('\t(90000.0 '+T_sources+')\n')
                        f.write(');\n\n')
                    print('Wrote out '+fpath)
            
                else:
                    params_copy['temperatureSourceType'] = 'given'
                    params_copy['temperatureForcing'] = \
"""sourceHeightsTemperature
(
    $windHeight
);

sourceTableTemperature
(
    (    0.0 0.0)
    (90000.0 0.0)
);"""

                tpl.copy_and_update(srcdir,dpath,
                                    os.path.join('constant','ABLProperties'),
                                    params_copy)


            # solver files
            os.makedirs(os.path.join(dpath,'system','sampling'))
            for fname in system_files:
                shutil.copy2(os.path.join(srcdir,'system',fname),
                             os.path.join(dpath,'system'))
            for fname in sampling_files:
                #shutil.copy2(os.path.join(srcdir,'system','sampling',fname),
                #             os.path.join(dpath,'system','sampling'))
                tpl.copy_and_update(srcdir,dpath,
                                    os.path.join('system','sampling',fname),
                                    params_copy)
            shutil.copy2(os.path.join(srcdir,'system','changeDictionaryDict.updateBCs.'+self.wdirname),
                         os.path.join(dpath,'system'))
            tpl.copy_and_update(srcdir,dpath,
                                os.path.join('system','setFieldsABLDict'),
                                params_copy)

            # scripts
            tpl.copy_and_update(srcdir,dpath,
                                'runscript.preprocess',
                                params_copy)

            tpl.copy_and_update(srcdir,dpath,
                                'runscript.solve.1',
                                params_copy)

        print('Done generating case directory in '+dpath)


    #--------------------------------------------------------------------------
    # sanity check routines
    #--------------------------------------------------------------------------

    def _alert(self,msg):
        messagebox.showwarning('Sanity check!',msg)

    def check_sanity(self):
        """TODO: add sanity checks go here"""

        # Check decomposition
        if (self.avgCellsPerCore > 80000):
            self._alert(str(self.avgCellsPerCore)+' cells/core... simulation will be slow')
        elif (self.avgCellsPerCore < 20000):
            self._alert(str(self.avgCellsPerCore)+' cells/core... maybe you need less cores')

        if self.decompTypeVar.get() == 'simple':
            decomp = self._text_to_list(self.decompOrder.get(), dtype=int)
            decompcores = int(np.prod(decomp))
            if not decompcores == int(self.nCores.get()):
                self._alert('"Simple" decomposition for '+str(decompcores)+
                            'cores, but '+self.nCores.get()+' cores expected')

        # Check domain
        if not (self.dx == self.dy == self.dz):
            self._alert('Did you mean to specify different spacings in each direction?')

        # Check surface conditions
        if (float(self.z0.get()) <= 0.0):
            self._alert('Surface roughness z0 <= 0 is likely to cause problems'+
                        ' for initialization, Rwall, and possibly qwall.')

        # Check IC/BCs
        if (self.velocityInitTypeVar.get() == 'table') \
                and (not self.sourceTypeVar.get() == 'profile'):
            self._alert('Need specified profile table for initialization')

        if (self.surfaceBCTypeVar.get() == 'fixed heating rate') \
                and (self.heating > 0):
            # heat flux INTO domain, with specified dT/dt
            self._alert('Fixed heating rate BC specified with heating, consider using fixed flux')
        elif (self.surfaceBCTypeVar.get() == 'fixed flux') \
                and (self.flux > 0):
            # heat flux OUT of domain (cooling)
            self._alert('Fixed flux BC specified with cooling, consider using fixed heating rate')
        elif (self.surfaceBCTypeVar.get() == 'fixed flux') \
                and (self.Tgrad_strong < 5):
            # heat flux INTO domain, with fixed flux but weak inversion
            self._alert('Convective conditions without a strong inversion layer,'+
                        ' {} K/(100 m) specified'.format(self.Tgrad_strong))

        if (self.temperatureInitTypeVar.get() == 'table'):
            listlist = self._text_to_listlist(self.profileTable.get('1.0',tk.END))
            data = np.array(listlist)
            if len(data) == 0:
                print('No specified profile')
            else:
                z = data[:,0]
                T = data[:,3]
                Tgrad = float(self.TGradUpper.get())
                dTdz_top = (T[-1]-T[-2]) / (z[-1]-z[-2])
                print('specified, calculated TGradUpper:',Tgrad,dTdz_top)
                #if not (Tgrad == dTdz_top):
                if not 1000*np.abs(Tgrad - dTdz_top) < 1e-8:
                    self._alert('Inconsistent upper boundary condition for temperature:\n' +
                            '\t{:.1f} K/km specified on boundary\n'.format(1000*Tgrad) + 
                            '\t{:.1f} K/km from specified interior profile'.format(1000*dTdz_top))

        # Background driving conditions
        if self.sourceTypeVar.get() == 'profile':
            wdir_variation = np.max(self.WD) - np.min(self.WD)
            if wdir_variation > 180.0:
                wdir_variation = np.min(self.WD) + 360. - np.max(self.WD)
            if wdir_variation > 90.0:
                self._alert('Veer results in a total wind direction change of'+
                            ' {} deg, need to use generalized inlet-outlet BC'.format(wdir_variation))


    #--------------------------------------------------------------------------
    # visualization routines
    #--------------------------------------------------------------------------

    # Note: importing matplotlib at the beginning causes a crash

    def plot_initial_conditions(self):
        self.update_profiles()
        inittype = self.velocityInitTypeVar.get()
        if inittype == 'geostrophic':
            U0 = float(self.U0Mag.get()) * np.ones(self.z_U.shape)
        elif inittype == 'log':
            zinv = float(self.zInversion.get())
            z0 = float(self.z0.get())
            Ug = float(self.U0Mag.get())
            kappa = self.params['kappa']
            ustar = kappa * Ug / np.log(zinv/z0)
            print('u* = {:g} m/s'.format(ustar))
            U0 = ustar / kappa * np.log(self.z_U/z0)
            U0[self.z_U > zinv] = Ug
        elif inittype == 'table':
            U0 = np.sqrt(self.U**2 + self.V**2)

        if DEBUG: print('Plotting initial conditions')
        import matplotlib.pyplot as plt
        plt.style.use(pltstyle)
        fig,ax = plt.subplots(ncols=2,sharey=True)
        fig.suptitle('Initial Conditions')
        ax[0].plot(U0, self.z_U)
        ax[1].plot(self.T, self.z_T)
        if not inittype == 'table':
            wdir = float(self.dir.get())
            ax[0].text(0, 1, r'direction: {:.1f}$^\circ$'.format(wdir),
                       horizontalalignment='left',
                       verticalalignment='top',
                       transform=ax[0].transAxes)
        ax[0].set_ylabel('height [m]')
        ax[0].set_xlabel('velocity [m/s]')
        ax[1].set_xlabel('temperature [K]')
        plt.show()

    def plot_background_conditions(self):
        self.update_profiles()
        listlist = self._text_to_listlist(self.profileTable.get('1.0',tk.END))
        data = np.array(listlist)
        z = data[:,0]
        U = data[:,1]
        V = data[:,2]
        T = data[:,3]

        wspd = np.sqrt(U*U + V*V)
        wdir = 180./np.pi * np.arctan2(-U,-V)
        wdir[wdir < 0] += 360.

        velocityInitType = self.velocityInitTypeVar.get()
        temperatureInitType = self.temperatureInitTypeVar.get()
        sourceType = self.sourceTypeVar.get()
        windtitle = ''
        temptitle = ''
        if velocityInitType == 'table':
            windtitle += 'initial conditions\n'
        if temperatureInitType == 'table':
            temptitle += 'initial conditions\n'
        if sourceType == 'profile':
            windtitle += 'background conditions'
            temptitle += 'background conditions'
        windtitle = windtitle.strip()
        temptitle = temptitle.strip()

        if DEBUG: print('Plotting profile table')
        import matplotlib.pyplot as plt
        plt.style.use(pltstyle)
        fig,ax = plt.subplots(ncols=3,sharey=True)
        fig.suptitle('Profile Table')
        ax[0].plot(wspd, z)
        ax[1].plot(wdir, z)
        ax[2].plot(T, z)
        ax[0].set_ylabel('height [m]')
        ax[0].set_xlabel('wind speed [m/s]')
        ax[1].set_xlabel('wind direction [deg]')
        ax[2].set_xlabel('temperature [K]')
        ax[0].set_title(windtitle,fontsize='small')
        ax[1].set_title(windtitle,fontsize='small')
        ax[2].set_title(temptitle,fontsize='small')
        overlaytext = 'NOT USED'
        if not sourceType == 'profile':
            if not velocityInitType == 'table':
                ax[0].text(0.5, 0.5, overlaytext, withdash=True,
                           horizontalalignment='center',
                           verticalalignment='center',
                           transform=ax[0].transAxes)
                ax[1].text(0.5, 0.5, overlaytext, withdash=True,
                           horizontalalignment='center',
                           verticalalignment='center',
                           transform=ax[1].transAxes)
            if not temperatureInitType == 'table':
                ax[2].text(0.5, 0.5, overlaytext, withdash=True,
                           horizontalalignment='center',
                           verticalalignment='center',
                           transform=ax[2].transAxes)
        plt.show()


#==============================================================================
#
# Code execution begins here
#

root = tk.Tk()  # the root window
#root.geometry('800x600')

# create instance of GUI
my_gui = MainWindow(root)
root.title('SOWFA Precursor Setup')

root.mainloop()
