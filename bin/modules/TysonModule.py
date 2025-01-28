# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 20:14:04 2025

@author: Arnab
"""
import numpy as np



def system_of_odes(t, y, params):
    dydt = np.zeros(6)
    
    # params_new = []
    
    # for param in params:
    #     param_new = param + 0.001*np.random.normal(0,np.sqrt(param))
    #     params_new.append(param_new)
    
    # params = np.array(params_new)
    
    Y = y[0]
    YP = y[1]
    C2 = y[2]
    CP = y[3]
    M = y[4]
    pM = y[5]
    
    k1norm = params[0]
    k2 = params[1]
    k3 = params[2]
    k4prime = params[3]
    k4 = params[4]
    k5 = params[5]
    k6 = params[6]
    k7 = params[7]
    k8 = params[8]
    k9 = params[9]
    
    k1 = k1norm*(C2 + CP + pM + M)
    CT = C2 + CP + pM + M
    
    
    dY = k1 - k2*Y - k3*CP*Y
    dYP = -k7*YP + k6*M
    dC2 = k6*M + k9*CP - k8*C2
    dCP = k8*C2 - k9*CP - k3*CP*Y
    dM = -k6*M -k5*M + pM*(k4prime + k4*(M/CT)**2)
    dpM = k3*CP*Y - pM*(k4prime + k4*(M/CT)**2) + k5*M
    
    
    dydt[0] = dY
    dydt[1] = dYP
    dydt[2] = dC2
    dydt[3] = dCP
    dydt[4] = dM
    dydt[5] = dpM
    return dydt