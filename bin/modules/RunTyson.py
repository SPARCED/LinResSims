# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 20:12:08 2025

@author: Arnab
"""

import numpy as np
from scipy.integrate import solve_ivp
from modules.TysonModule import system_of_odes


def RunTyson(y0,th,params):
    
    t_min = th*60
    
    t_span = (0,t_min)
    
    n_tp = int(t_span[1]*2+1)
    
    solution = solve_ivp(system_of_odes, t_span, y0, t_eval=np.linspace(t_span[0], t_span[1], n_tp), args=(params,),method='LSODA')
    t = solution.t
    y = solution.y
    
    tout_all = t
    xoutS_all = np.transpose(y)

    
    return xoutS_all, tout_all