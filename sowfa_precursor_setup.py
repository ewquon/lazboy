#!/usr/bin/env python
#
# GUI for semi-automatic SOWFA precursor set up
#
# Written by Eliot Quon (eliot.quon@nrel.gov), 2018-07-09
#
from __future__ import print_function
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText

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

        self.init_vars()
        self.init_window()
        self.restore_defaults()


    def init_vars(self):
        # Need to initialize some variables before we can create the widgets
        # However, most values associated with Entry boxes will be created
        # on the fly.
        self.lastrow = -1

        """decomposition controls"""
        #self.nCores = 216  # Number of cores on which to run this case.
        #self.decompType = 'simple'  # Decomposition algorithm.  "simple" and "scotch" are good choices.
        self.decompTypeVar = tk.StringVar(self.master, value='simple')
        #self.decompOrder = [6,6,6]  # Order of the decomposition number of partitions in (x y z)-directions.

        """general conditions"""
        #self.TRef = 300.0  # Reference potential temperature (K).
        self.coriolisVar = tk.IntVar(self.master, value=1)
        #self.latitude = 40.0  # Latitude on the Earth of the site (deg).
        #self.EarthPeriod = 24.0  # Earth's rotation period (hr).

        """atmosphere controls"""
        self.velocityInitTypeVar = tk.StringVar(self.master, value='geostrophic')
        self.temperatureInitTypeVar = tk.StringVar(self.master, value='simple')
#        self.U0Mag = 8.0  # Initial condition for wind speed (m/s).
#        self.dir = 270.0  # Initial condition for wind direction (deg).
#        self.windHeight = 80.0  # Height at which to drive mean wind to U0Mag/dir (m).
#        self.p_rgh0 = 0.0  # Initial pressure (minus the hydrostatic variation and normalized by density) (m^2/s^2).
#        self.nuSgs0 = 0.0  # Initial SGS viscosity (m^2/s).
#        self.k0 = 0.1  # Initial SGS turbulent kinetic energy (m^2/s^2).
#        self.kappat0 = 0.0  # Initial SGS temperature diffusivity (m^2/s).
#        self.TGradUpper = 0.003  # Potential temperature gradient above the strong inversion (K/m).
#        self.zInversion = 750.0  # Height of the middle of the initial strong capping inversion (m).
#        self.inversionWidth = 100.0  # Vertical width of the intial strong capping inversion (m).
#        self.TBottom = 300.0  # Initial potential temperature at bottom of strong capping inversion (K).
#        self.TTop = 305.0  # Initial potential temperature at top of strong capping inversion (K).
#
#        """surface controls"""
#        self.qwall = [0.0,0.0,0.0]  # Temperature flux at wall (modify the z-value).  A negative value is flux into domain (K-m/s).
#        self.Rwall = [0.0,0.0,0.0,0.0,0.0,0.0]  # Initial wall shear stress (m^2/s^2).
#        self.kappa = 0.4  # von Karman constant.
#        self.z0 = 0.15  # Surface roughness (m).
#        self.heatingRate = 0.0  # Surface temperature change rate (when not directly setting temperature flux) (K/s).
#
#        """advanced controls"""
#        # surface conditions
#        self.wallModelAverageType = 'planarAverage'  # Treat surface stress wall model locally ("local") or with planar averaging ("planarAverage").
#        self.betaM = 16.0  # Monin-Obukhov wall shear stress model constant.
#        self.gammaM = 5.0  # Monin-Obukhov wall shear stress model constant.
#        self.betaH = 9.0  # Monin-Obukhov wall temperature flux model constant.
#        self.gammaH = 7.8  # Monin-Obukhov wall temperature flux model constant.
#        self.alphaH = 1.0  # Monin-Obukhov wall temperature flux model constant.
#
#        # planar averaging and source term statistics options.
#        self.statisticsOn = True  # Gather planar-averaged flow statistics.
#        self.statisticsFrequency = 5  # Frequency in time steps of statistics gathering.
#
#        # transport properties
#        self.Pr = 0.7  # Molecular Prandtl number.
#        self.Prt = 0.33333333  # Turbulent Prandtl number.
#        self.nu = 1.0E-5  # Molecular viscosity (m^2/s).
#
#        # SGS model inputs
#        self.LESModel = 'oneEqEddyABL'  # SGS model selection.
#        self.ce = 0.93  # SGS model constant.
#        self.ck = 0.0673  # SGS model constant.

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
        # top
        text = tk.Label(self.master, text='Canonical ABL template:')
        #stability_list = ['unstable','neutral','stable']
        template_list = tpl.get_templates()
        self.template_option = tk.OptionMenu(self.master,
                                             self.template,
                                             *template_list)
        self.restore_button = tk.Button(self.master,
                                        text='Restore defaults',
                                        command=self.restore_defaults)
        self.save_button = tk.Button(self.master,
                                     text='Save template',
                                     command=self.save_template)
        text.grid(row=self.nextrow(), sticky=tk.E)
        self.template_option.grid(row=self.lastrow, column=1, sticky=tk.W)
        self.restore_button.grid(row=self.lastrow, column=2)
        self.save_button.grid(row=self.lastrow, column=3)

        # main part of window
        self.init_domain_controls()
        self.init_decomp_controls()
        self.init_general_controls()
        self.init_atmosphere_controls()
        self.init_surface_controls()
        self.init_advanced_controls()

        # bottom
        self.generate_button = tk.Button(self.master,
                                         text='Generate case files!',
                                         command=self.generate)
        self.generate_button.grid(row=self.nextrow(), columnspan=4,
                                  sticky='nesw') # fill span


    def init_domain_controls(self):
        section = tk.LabelFrame(self.master, text='Domain')
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
        section = tk.LabelFrame(self.master, text='Domain Decomposition')
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
        section = tk.LabelFrame(self.master, text='General Settings')
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
        section = tk.LabelFrame(self.master, text='Atmospheric Conditions')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)

        header = tk.Label(section, text='Initial conditions:')
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
        self.U0Mag = self.EntryRow(section,
                                   'U0Mag',
                                   'Initial condition for wind speed (m/s).')
        self.dir = self.EntryRow(section,
                                 'dir',
                                 'Initial condition for wind direction (deg).')
        self.windHeight = self.EntryRow(section,
                                        'windHeight',
                                        'Height at which to drive mean wind to U0Mag/dir (m).')
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
                                        'Height of the middle of the initial strong capping inversion (m).')
        self.inversionWidth = self.EntryRow(section,
                                            'inversionWidth',
                                            'Vertical width of the intial strong capping inversion (m).')
        self.TBottom = self.EntryRow(section,
                                     'TBottom',
                                     'Initial potential temperature at bottom of strong capping inversion (K).')
        self.TTop = self.EntryRow(section,
                                  'TTop',
                                  'Initial potential temperature at top of strong capping inversion (K).')

        header = tk.Label(section, text='Background conditions:')
        header.grid(row=self.nextrow(), column=0, sticky=tk.W)
        self.plot_bkgnd_button = tk.Button(section,
                                           text='PLOT',
                                           command=self.plot_background_conditions)
        self.plot_bkgnd_button.grid(row=self.lastrow, column=3, sticky=tk.E)

        """Optional profile table for finer control over the background
        and initial conditions. Note, I'm assuming that if a profile is
        specified for the background conditions, the same profile should
        be used to initialize the solution with setFieldsABL.
        """
        namestr = tk.Label(section, text='profileTable')
        namestr.grid(row=self.nextrow(), column=0)
        text = ScrolledText(section)
        text.grid(row=self.lastrow, column=1, columnspan=2, sticky='ew')
        self.profileTable = text


    def init_surface_controls(self):
        section = tk.LabelFrame(self.master, text='Surface Conditions')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)


    def init_advanced_controls(self):
        section = tk.LabelFrame(self.master, text='Advanced Controls')
        section.grid(row=self.nextrow(), columnspan=4, pady=5, sticky=tk.W)

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

        for name, val in params.iteritems():
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
        self.update_decomp(self.decompTypeVar.get())
        self.update_coriolis()
        self.update_velocity_init(self.velocityInitTypeVar.get())
        self.update_temperature_init(self.temperatureInitTypeVar.get())


    def save_template(self):
        print('Saving template')


    def generate(self):
        print('Generating case files')


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
            self.dxText['text'] = 'dx = {} m'.format((x1-x0)/nx)
            self.dyText['text'] = 'dy = {} m'.format((y1-y0)/ny)
            self.dzText['text'] = 'dz = {} m'.format((z1-z0)/nz)
            self.NcellsText['text'] = 'Ncells = {:g}'.format(nx*ny*nz)

            self.coresNeeded = int(self.nCores.get())/int(self.PPN.get())
            self.NnodesText['text'] = '{} compute nodes'.format(self.coresNeeded)

            avgCellsPerCore = int(nx*ny*nz/float(self.nCores.get()))
            self.cellsPerCoreText['text'] = 'average {} cells/core'.format(avgCellsPerCore)

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

    def update_velocity_init(self,initType):
        print('velocityInitType = '+initType)

    def update_temperature_init(self,initType):
        print('temperatureInitType = '+initType)


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

def on_configure(event):
    # update scrollregion after starting 'mainloop'
    # when all widgets are in canvas
    canvas.configure(scrollregion=canvas.bbox('all'))

#------------------------------------------------------------------------------

root = tk.Tk()  # the root window
#root.wm_geometry(...)

# create canvas with scrollbar
canvas = tk.Canvas(root)
canvas.pack(side=tk.LEFT)
scrollbar = tk.Scrollbar(root, command=canvas.yview)
scrollbar.pack(side=tk.LEFT, fill='y')

canvas.configure(yscrollcommand=scrollbar.set)

# update scrolling region after starting main loop
canvas.bind('<Configure>', on_configure)

# put frame in canvas
#frame = tk.Frame(root, relief=tk.GROOVE)
#frame = tk.Frame(canvas)

# create instance of window
#my_gui = MainWindow(root)

my_gui = MainWindow(canvas)
root.title('SOWFA Precursor Setup')

canvas.create_window((0,0), window=my_gui, anchor='nw')

root.mainloop()
