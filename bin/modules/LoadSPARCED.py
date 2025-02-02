# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 18:37:48 2025

@author: Arnab
"""

import numpy as np
import sys
import importlib
import os


def LoadSPARCED(sim_config,wd):
    
    sbml_file = "SPARCED.xml"
    model_name= sbml_file[0:-4]
    model_output_dir = model_name
    sys.path.insert(0, os.path.join(wd,model_output_dir))
    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver()
    solver.setMaxSteps = 1e10
    ts = 30

    species_all = list(model.getStateIds())
    model.setTimepoints(np.linspace(0,ts))
    STIMligs_id = ['E', 'H', 'HGF', 'P', 'F', 'I', 'INS']
    species_initializations = np.array(model_module.getModel().getInitialStates())
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
    
    for l,lig in enumerate(STIMligs_id):
        species_initializations[species_all.index(lig)] = 0

    kwargs_default = {}
    kwargs_default['flagD'] = 0
    kwargs_default['th'] = sim_config["preinc_time"]
    kwargs_default['spdata'] = species_initializations
    kwargs_default['genedata'] = []
    kwargs_default['sbml_file'] = sbml_file
    kwargs_default['model'] = model
    
    cc_marker = sim_config['cc_marker']
    
    model_specs = {'species_all': species_all, 'cc_marker': cc_marker}
    return model_specs, kwargs_default