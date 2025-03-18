# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 20:12:08 2025

@author: Arnab
"""

import numpy as np
from scipy.integrate import solve_ivp
from modules.TysonModule import system_of_odes
import math

def round_tp(num):
    
    remainder = num - math.floor(num)
    if remainder > 0 and remainder < 0.5:
        remainder = 0
    elif remainder >= 0.5 and remainder < 1.0:
        remainder = 0.5
    
    new_num = math.floor(num) + remainder
    
    return new_num


def RunTyson(th,spdata,params):
    
    t_min = th*60
    
    t_span = (0,t_min)
    
    n_tp = int(t_span[1]*2+1)
    
    tp_cell = np.linspace(t_span[0], t_span[1], n_tp)
    
    tp_rounded = [round_tp(tp) for tp in tp_cell]
    
    solution = solve_ivp(system_of_odes, t_span, spdata, t_eval=tp_rounded, args=(params,),method='LSODA')
    t = solution.t
    y = solution.y
    
    tout_all = t
    xoutS_all = np.transpose(y)

    
    return xoutS_all, tout_all