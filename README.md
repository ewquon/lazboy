# lazboy
A cushy set of GUI-based tools to help with setting up and running SOWFA
(Simulator fOr Wind Farm Applications: https://github.com/NREL/SOWFA)
simulations. 

These scripts depend on the Tkinter package for the GUI and are written with
both Python 2 and 3 compatibility in mind. 

The scripts should be from the command line in a terminal with X11 forwarding
(either locally or remotely), or from a virtual desktop environment (e.g.,
FastX on Peregrine: https://www.nrel.gov/hpc/software-fastx.html):
```
python /path/to/lazboy/sowfa_precursor_setup.py
```
If the GUI is laggy and you're running remotely, then FastX is probably your
best option.

## Getting started: sowfa_precursor_setup
* FIRST RUN: You will be prompted for your email and primary HPC allocation
(for generating runscripts with email notification). This will be stored in
`$HOME/.sowfa_defaults`. Depending on your screen configuration, the prompts
may be hidden behind the main window. 
* The default canonical precursor template is "neutral". Available precursor
templates are provided in the dropdown menu at the top of the window, and
populated from `/path/to/lazboy/simulation_templates/`. These are based on Matt
Churchfield's setUp file and saved in yaml format. You can add manually create
new templates or save the current configuration as a template using the "Save
template" button at the top of the screen.
* In addition to the default template, this utility will also look in your
current working directory for additional yaml files.
* When you're done configuring your simulation, click "Generate case files!" at
the bottom of the screen. This will prompt you for a directory into which to
create the simulation files ("constant", "system", "0.original", "setUp", and
"runscript.\*"). A yaml file describing the configuration will also be saved
into `/path/to/new/simulation/setup.yaml`; this can directly be used as a
template. 
* Sanity checks are performed on-the-fly when you interact with some of the
  configuration options, and also when you click "Generate case files!". 

![SOWFA precursor setup screenshot](https://raw.githubusercontent.com/ewquon/lazboy/master/screenshot.png)

### quirks
* In linux, the GUI rendering looks super old-school. However the functionality
should be mostly the same for Mac/Windows/Linux. 
* Scrolling only works if you move your mouse over the scrollbar. 
* In linux, the "Specify new simulation directory" box doesn't have a make-new-
directory button, so instead enter the name you want and click "OK" twice.

If something is behaving unexpectedly, check the terminal for some (hopefully)
useful verbose degugging information...

