# -*- coding: utf-8 -*-
"""
Created on Mon May 13 12:13:41 2024

@author: Arnab
"""

# import required libraries

import pickle
import re

import itertools
import math

import os
import sys

import libsbml
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
import copy

from Bio import Phylo

from io import StringIO

from scipy.interpolate import interp1d
from scipy.stats import percentileofscore
from scipy.stats import gaussian_kde
import math
import seaborn as sns
import itertools
import plotly.figure_factory as ff
import plotly.io as pio


import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.sans-serif'] = ['Arial']


#%%
cd = os.getcwd()
wd = os.path.dirname(cd)
sys.path.append(os.path.join(wd,'bin'))

sbml_file = "SPARCED.xml"


sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(os.path.join(wd,sbml_file))
sbml_model = sbml_doc.getModel()

species_all = [str(x.getId()) for x in list(sbml_model.getListOfSpecies())]





#%%
output_dir_main = os.path.join(wd,'output')

exp_title = 'in_silico_drs'
output_main = os.path.join(wd,'output',exp_title)



if not os.path.exists(os.path.join(wd, 'output', 'in_silico_drs_summary')):
    os.mkdir(os.path.join(wd, 'output', 'in_silico_drs_summary'))


dir_doses_all = os.listdir(os.path.join(output_main,'drs_alpel','drs_alpel_rep1'))

doses_all = [float(x.split('_')[-1]) for x in dir_doses_all]

doses_all.sort()


#%%

# import class for reading dose response outputs

from modules.drsPlotting import drs_dict


#

#%% dose response population dynamics



def drs_summarize(drug,dose_lvl,output_dir=output_main): # measures alive cell count from each simulation
    
    drs_dose = {}
    for rep in range(10):
        print('now running...'+str(drug)+'...'+str(dose_lvl)+'...'+str(rep+1))
        drs_dict0 = drs_dict(output_dir,drug,rep+1,dose_lvl)
        pd,tp,td = drs_dict0.pop_dyn()
        drs_rep = {}
        drs_rep['cellpop'] = pd
        drs_rep['tout'] = tp
        drs_rep['t_death'] = td
        drs_dose['r'+str(rep+1)] = drs_rep
    return drs_dose



#%%



def drs_death_summarize(drug_exp,dose_lvl,output_dir=output_main,doses_all=doses_all): # measures alive cell count from each simulation
    
    drug = drug_exp[:5].lower()
    
    df_death_all = pd.DataFrame(columns=['drug','dose','rep','t_death'])


    for rep in range(10):
        print('now running...'+str(drug)+'...'+str(dose_lvl)+'...'+str(rep+1))
        drs_dict0 = drs_dict(output_dir,drug,rep+1,dose_lvl)
        popd,tp,td = drs_dict0.pop_dyn()

        
        td = td[~np.isnan(td)]
        
        df_death = pd.DataFrame(columns=['drug','dose','rep','t_death'])
        df_death['drug'] = [drug_exp]*len(td)
        df_death['dose'] = [doses_all[dose_lvl]]*len(td)
        df_death['rep'] = ['rep'+str(rep+1)]*len(td)
        df_death['t_death'] = td
        
        df_death_all = pd.concat((df_death_all,df_death),axis=0)
        
        
    return df_death_all







#%% generate population dynamics summary datasets



drugs_exp = ['Alpelisib','Neratinib','Trametinib','Palbociclib']



#%

for dr_idx in range(len(drugs_exp)):
    
    drug = drugs_exp[dr_idx]
    
    df0 = pd.DataFrame(columns=['drug','dose','rep','t_death'])
    
    drs_drug = {}
    
    for d in range(10):
        
        df_death_dose = drs_death_summarize(drug, d)
        df0 = pd.concat((df0,df_death_dose),axis=0)
        
    df0.to_csv(os.path.join(wd,'output','in_silico_drs_summary',f'drs_death_{drug[:5].lower()}.txt'),sep='\t',index=False,header=True)   


