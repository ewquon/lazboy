from __future__ import print_function
import os
import yaml

mypath = os.path.dirname(os.path.realpath(__file__))

def get_templates():
    dpath = os.path.join(mypath, 'simulation_templates')
    listing = os.listdir(dpath)
    templates = [ os.path.splitext(name)[0] for name in listing
                    if os.path.isfile(os.path.join(dpath,name))
                        and name.endswith('.yaml') ]
    return templates

def read_template(name):
    fpath = os.path.join(mypath, 'simulation_templates', name+'.yaml')
    print('Loaded template: '+fpath)
    with open(fpath,'r') as f:
        params = yaml.load(f)
    return params

def save_template(name):
    fpath = os.path.join(mypath, 'simulation_tempolates', name)
    if os.path.isfile(fpath):
        print('Template '+name+' already exists! No file was written.')
        return
