# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 19:55:09 2025

@author: Arnab
"""
import numpy as np
import re


def assign_tasks(rank,n_cells,size):
    
    cells_per_rank = n_cells // int(size)
    remainder = n_cells % int(size)
    
    if rank < remainder:
        my_cells = cells_per_rank + 1
        start_cell = rank * my_cells + 1
    else:
        my_cells = cells_per_rank
        start_cell = rank * cells_per_rank + remainder + 1
        
    return start_cell, start_cell + my_cells


def args_override(config,args): # override sim_config options with commandline arguments
    config_new = config
    
    args_dict = vars(args)  # Convert argparse Namespace to dictionary
    for key, value in args_dict.items():
        if value is not None and key != "sim_config":  # Skip None values and 'config' key
            if key in config_new.keys():  # Top-level keys
                # # Handle boolean conversion for enable_logging if necessary
                # config[key] = value.lower() == "true" if key == "enable_logging" else value
                config_new[key] = value
            elif key in config.get("stim", {}):  # Nested keys in 'parameters'
                config["stim"][key] = value
    
    
    # if args.cellpop is not None:
    #     config_new['cellpop'] = int(args.cellpop)
    # if args.exp_time is not None:
    #     config_new['exp_time'] = float(args.exp_time)

    return config_new


def evaluate_formula(values, labels, formula, locals_dict):
    """
    Evaluates an algebraic expression using values mapped by labels.
    
    Parameters:
    values (np.array): Array of values corresponding to labels.
    labels (list): List of labels corresponding to indices in `values`.
    formula (str): Algebraic expression as a string using label names.
    
    Returns:
    float: Evaluated result.
    """
        
    sp_formula = re.findall(r'[a-zA-Z]\w*',formula)
    sp_formula = list(np.unique(sp_formula))
    
    if 'e' in sp_formula:
        sp_formula.remove('e')
        

    result = None
        
    for i in range(len(sp_formula)):
        exec(f"{sp_formula[i]} = values[labels.index('{sp_formula[i]}')]",globals(), locals_dict)
        
    exec(f"flagA = {formula}",globals(), locals_dict)
    
    return result
