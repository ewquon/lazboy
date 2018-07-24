#!/usr/bin/env python
#
# GUI for semi-automatic SOWFA precursor set up
#
# Written by Eliot Quon (eliot.quon@nrel.gov), 2018-07-09
#
from __future__ import print_function
import os
import yaml
import numpy as np

import tkinter as tk
from tkinter.scrolledtext import ScrolledText
try:
    # Python 2.7
    import tkFileDialog as filedialog
    import tkSimpleDialog as simpledialog
except ImportError:
    # Python 3
    import tkinter.filedialog as filedialog
    import tkinter.simpledialog as simpledialog

import template as tpl

DEBUG = True

#class MainWindow(object):
class MainWindow(tk.Frame):
    """GUI for configuring the setUp file and generating a SOWFA
    precursor input deck
    """

    def __init__(self, master, default_template='neutral'):
        tk.Frame.__init__(self)

        self.master = master
        self.template = tk.StringVar(master, value=default_template)
        self.params = tpl.read_template(default_template)

        self.user_info()

        self.init_vars()
        self.init_window()
        self.restore_defaults()


    def user_info(self):
        """Read or set up (first time) system information"""
        configpath = os.path.join(os.environ['HOME'],'.sowfa_defaults')
        try:
            with open(configpath,'r') as f:
                config = yaml.load(f)
        except IOError:
            config = dict()
            config['email'] = simpledialog.askstring(
                                    title='Default SOWFA configuration',
                                    prompt='Notification email',
                                    initialvalue='Richard.Astley@nrel.gov') 
            config['allocation'] = simpledialog.askstring(
                                    title='Default SOWFA configuration',
                                    prompt='Default allocation',
                                    initialvalue='windsim') 
            with open(configpath,'w') as f:
                yaml.dump(config, f, default_flow_style=False)
            print('Saved default SOWFA configuration info to '+configpath)
        if DEBUG: print(config)

            
    def init_vars(self):
        # Need to initialize some variables before we can create the widgets
        # However, most values associated with Entry boxes will be created
        # on the fly.
        self.lastrow = -1

        """decomposition controls"""
        self.decompTypeVar = tk.StringVar(self.master, value='simple')

        """general conditions"""
        self.coriolisVar = tk.IntVar(self.master, value=1)

        """atmosphere controls"""
        self.velocityInitTypeVar = tk.StringVar(self.master, value='geostrophic')
        self.temperatureInitTypeVar = tk.StringVar(self.master, value='simple')
        self.sourceTypeVar = tk.StringVar(self.master, value='single height')
        #self.p_rgh0 = 0.0  # Initial pressure (minus the hydrostatic variation and normalized by density) (m^2/s^2).
        #self.nuSgs0 = 0.0  # Initial SGS viscosity (m^2/s).
        #self.k0 = 0.1  # Initial SGS turbulent kinetic energy (m^2/s^2).
        #self.kappat0 = 0.0  # Initial SGS temperature diffusivity (m^2/s).
 
        """surface controls"""
        self.surfaceBCTypeVar = tk.StringVar(self.master, value='fixed flux')
        #self.Rwall = [0.0,0.0,0.0,0.0,0.0,0.0]  # Initial wall shear stress (m^2/s^2).
        #self.kappa = 0.4  # von Karman constant.
 

    #==========================================================================
    #
    # TK wrappers
    #

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
        entry = tk.Checkbutton(master, text=description, variable=variable,
                               **kwargs)
        #descstr = tk.Label(master, text=description)
        namestr.grid(row=self.nextrow(), column=0)
        entry.grid(row=self.lastrow, column=1, sticky='ew')
        #descstr.grid(row=self.lastrow, column=2, sticky=tk.W)
        return entry

    #==========================================================================
    #
    # GUI setup
    #

    def init_window(self):
        self.top = tk.Canvas(self.master)
        self.middle = tk.Canvas(self.master)
        self.bottom = tk.Canvas(self.master)
        self.top.pack(side='top', fill='x')
        self.middle.pack(side='top', fill='both', expand=True)
        self.bottom.pack(side='top', fill='x')

        # top
        text = tk.Label(self.top, text='Canonical ABL template:')
        template_list = tpl.get_templates()
        self.template_option = tk.OptionMenu(self.top,
                                             self.template,
                                             *template_list)
        self.restore_button = tk.Button(self.top,
                                        text='Restore defaults',
                                        command=self.restore_defaults)
        self.save_button = tk.Button(self.top,
                                     text='Save template',
                                     command=self.save_template)
        #text.grid(row=0, sticky=tk.E)
        #self.template_option.grid(row=0, column=1, sticky=tk.W)
        #self.restore_button.grid(row=0, column=2)
        #self.save_button.grid(row=0, column=3)
        text.pack(side='left')
        self.template_option.pack(side='left', padx=20)
        self.restore_button.pack(side='left', padx=5)
        self.save_button.pack(side='left')

        # main part of window
        self.create_main_canvas()
        self.init_domain_controls()
        self.init_decomp_controls()
        self.init_general_controls()
        self.init_atmosphere_controls()
        self.init_surface_controls()
        self.init_advanced_controls()

        # bottom
        self.generate_button = tk.Button(self.bottom,
                                         text='Generate case files!',
                                         command=self.generate)
        #self.generate_button.grid(row=2), columnspan=4,
        #                          sticky='nesw') # fill span
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
                                        'Potential temperature gradient above the strong inversion (K/m).')
        self.zInversion = self.EntryRow(section,
                                        'zInversion',
                                        'Height of the middle of the initial strong capping inversion (m).',
                                        vcmd=vcmd)
        self.inversionWidth = self.EntryRow(section,
                                            'inversionWidth',
                                            'Vertical width of the intial strong capping inversion (m).')
        self.TBottom = self.EntryRow(section,
                                     'TBottom',
                                     'Initial potential temperature at bottom of strong capping inversion (K).')
        self.TTop = self.EntryRow(section,
                                  'TTop',
                                  'Initial potential temperature at top of strong capping inversion (K).')

        header = tk.Label(section,
                          text='Background conditions', font='-slant italic')
        header.grid(row=self.nextrow(), column=0, sticky=tk.W)
        self.plot_bkgnd_button = tk.Button(section,
                                           text='PLOT',
                                           command=self.plot_background_conditions)
        self.plot_bkgnd_button.grid(row=self.lastrow, column=3, sticky=tk.E)

        self.sourceType = self.OptionMenuRow(section, 'sourceType',
                                             'Computed source type',
                                             args=['single height','column'],
                                             variable=self.sourceTypeVar,
                                             command=self.update_source_type)
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
        namestr = tk.Label(section, text='profileTable')
        namestr.grid(row=self.nextrow(), column=0)
        text = ScrolledText(section)
        text.grid(row=self.lastrow, column=1, columnspan=2, sticky='ew')
        self.profileTable = text


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

    def init_advanced_controls(self):
        #section = tk.LabelFrame(self.master, text='Advanced Controls')
        section = tk.LabelFrame(self.scrollableframe, text='Advanced Controls',
                                font='-weight bold')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)

    #==========================================================================
    #
    # validation routines
    #

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
            dx = (x1-x0) / nx
            dy = (y1-y0) / ny
            dz = (z1-z0) / nz
            self.zsurf = z0
            self.dz = dz
            self.dxText['text'] = 'dx = {} m'.format(dx)
            self.dyText['text'] = 'dy = {} m'.format(dy)
            self.dzText['text'] = 'dz = {} m'.format(dz)
            self.NcellsText['text'] = 'Ncells = {:g}'.format(nx*ny*nz)

            self.coresNeeded = int(self.nCores.get())/int(self.PPN.get())
            self.NnodesText['text'] = '{} compute nodes'.format(self.coresNeeded)

            avgCellsPerCore = int(nx*ny*nz/float(self.nCores.get()))
            self.cellsPerCoreText['text'] = 'average {} cells/core'.format(avgCellsPerCore)


    def update_profiles(self,name=None):
        self.profileTable.delete('1.0', tk.END)
        self.profileTable.insert(tk.END, '//   z       U       V       T\n')
        nz = int(self.nz.get())
        self.z = self.zsurf + self.dz/2 + self.dz*np.arange(nz)
        zref = float(self.windHeight.get())
        Uref = float(self.U0Mag.get())
        wdir = 270. - float(self.dir.get())
        if wdir < 0: wdir += 360.
        wdir *= np.pi/180.
        alpha = float(self.alpha.get())
        veer = float(self.veer.get()) * np.pi/180.
        zinv = float(self.zInversion.get())
        self.WS = Uref * (self.z / zref)**alpha
        self.WD = wdir + veer/zref * (self.z - zref)
        freeatm = self.z >= zinv
        WDtop = self.WD[freeatm][0]
        self.WD[freeatm] = WDtop
        self.U = self.WS * np.cos(self.WD)
        self.V = self.WS * np.sin(self.WD)
        for zi,Ui,Vi in zip(self.z,self.U,self.V):
            rowdata = '{}  {}  {}\n'.format(zi,Ui,Vi)
            self.profileTable.insert(tk.END, rowdata)
        return True # so the widget isn't disabled

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


    def update_source_type(self,sourceType):
        if sourceType=='single height':
            self.profileTable.config(state='disabled',fg='light grey')
        elif sourceType=='column':
            self.profileTable.config(state='normal',fg='black')
        else:
            raise ValueError('Unexpected source type: '+sourceType)
        self.update_profiles()


    def update_surface_bc(self,name=None):
        R = 287. # [J/(kg-K)]
        pref = 101300. # standard atmospheric pressure at sea level [Pa]
        rhoref = pref / R / float(self.TRef.get())
        Cp = 3.5*R # gamma/(gamma-1) * R
        qwall = float(self.qwall.get())
        heatingRate = float(self.heatingRate.get())

        flux = qwall * rhoref * Cp
        fluxdesc = ''
        if flux < 0:
            fluxdesc = 'HEATING'
        elif flux > 0:
            fluxdesc = 'COOLING'
        self.fluxInfoText['text'] = '= {} W/m^2 {}'.format(np.abs(flux),fluxdesc)

        heatingdesc = ''
        if heatingRate > 0:
            heatingdesc = 'HEATING'
        elif heatingRate < 0:
            heatingdesc = 'COOLING'
        self.heatingInfoText['text'] = heatingdesc

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
        print('velocityInitType = '+initType)

    def update_temperature_init(self,initType):
        print('temperatureInitType = '+initType)


    #==========================================================================
    #
    # actions
    #

    def restore_defaults(self):
        """It's kind of a PITA, but binding tk.StringVar to the widgets
        and trying to udpate it in the validate command breaks the
        callback routine. So instead of being able to just call 
            tk.StringVar.set(value),
        we instead have to call
            tk.Entry.delete(0, END)
            tk.Entry.insert(0, value).
        """
        print('Restoring defaults for '+self.template.get())
        params = tpl.read_template(self.template.get())
        if DEBUG: print(params)
        self.params = params

        for name, val in params.items():
            try:
                widget = getattr(self, name)
            except AttributeError:
                print('Note: widget "'+name+'" does not exist')
            else:
                if widget.winfo_class() == 'Entry':
                    if DEBUG: print('Setting entry "'+name+'"')
                    widget.delete(0, tk.END)
                    widget.insert(0, val)
                else:
                    # we have a control variable instead of widget
                    if DEBUG: print('Setting variable for "{}" ({})'.format(name,widget.winfo_class()))
                    wvar = getattr(self, name+'Var')
                    wvar.set(val)
        
        """update active text widgets"""
        self.calc_grid_res()
        #self.update_profiles()
        self.update_source_type(self.sourceTypeVar.get()) # calls update_profiles()
        self.update_decomp(self.decompTypeVar.get())
        self.update_coriolis()
        self.update_velocity_init(self.velocityInitTypeVar.get())
        self.update_temperature_init(self.temperatureInitTypeVar.get())
        self.update_surface_bc()


    def get_all_params(self):
        for name, val in self.params.iteritems():
            dtype = type(val)
            try:
                widget = getattr(self, name)
            except AttributeError:
                pass
            else:
                if widget.winfo_class() == 'Entry':
                    if DEBUG: print('Updated variable "{}" ({}) from Entry'.format(name,str(dtype)))
                    wvar = widget.get()
                else:
                    # we have a control variable instead of widget
                    if DEBUG:
                        print('Updated variable "{}" ({}) from {}'.format(
                            name, str(dtype), widget.winfo_class()) )
                    wvar = getattr(self, name+'Var').get()

                if dtype is list:
                    dtype = type(val[0])
                    vals = wvar.strip('[]').split(',')
                    vals = [ dtype(val) for val in vals ]
                    self.params[name] = vals
                    if DEBUG: print('  reconstructed list: {}'.format(vals))
                else:
                    self.params[name] = dtype(wvar)


    def save_template(self,fpath=None):
        if fpath is None:
            self.get_all_params()
            fpath = filedialog.asksaveasfilename(
                        initialdir=os.path.join(tpl.mypath, 'simulation_templates'),
                        title='Save precursor configuration',
                        filetypes=(('yaml files','*.yaml'),('all files','*,*')),
                    )

        if not fpath == '':
            if DEBUG: print(tpl.yaml_template.format(**self.params))
            with open(fpath,'w') as f:
                f.write(tpl.yaml_template.format(**self.params))
            print('Wrote out '+fpath)

    def generate(self):
        self.get_all_params()
        self.check_sanity()

        dpath = filedialog.askdirectory(
                    initialdir=os.getcwd(),
                    title='Specify new simulation directory',
                )

        if not dpath == '':
            print('Generating case files in '+dpath)
            self.save_template(fpath=os.path.join(dpath,'config.yaml'))


    #==========================================================================
    #
    # sanity check routines
    #

    def check_sanity(self):
        print('TODO: add sanity checks!')

    #==========================================================================
    #
    # visualization routines
    #
    # Note: importing matplotlib at the beginning causes a crash

    def plot_initial_conditions(self):
        # TODO: This is a stub.
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(ncols=2)
        fig.suptitle('Initial Conditions')
        plt.show()

    def plot_background_conditions(self):
        # TODO: This is a stub.
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(ncols=2)
        fig.suptitle('Background Conditions')
        plt.show()


#------------------------------------------------------------------------------
# set up scrolling

#def on_configure(event):
#    # update scrollregion after starting 'mainloop'
#    # when all widgets are in canvas
#    print('update scrollregion')
#    canvas.configure(scrollregion=canvas.bbox('all'))

#------------------------------------------------------------------------------

root = tk.Tk()  # the root window
#root.geometry('800x600')

# make it cover the entire screen
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
#root.overrideredirect(1) # get rid of titlebar / menus
root.geometry('{:d}x{:d}+0+0'.format(w,h))

# create instance of window
my_gui = MainWindow(root)
root.title('SOWFA Precursor Setup')

root.mainloop()
