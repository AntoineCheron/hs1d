# -*- coding: utf-8 -*-
"""
Created on Tue May 23 10:27:08 2017

@author: Quentin Courtois
"""
import numpy as np
from Initial import Computation

def compute_Smax(f, w, soil_depth):
    Smax = f*np.reshape(w,(len(w),1)) * np.reshape(soil_depth,(len(soil_depth),1))
    return Smax

def compute_sin(percentage_loaded, Smax, Init, BC):
    if isinstance(Init, int):
        percentage_loaded = percentage_loaded
        if isinstance(percentage_loaded, float) or isinstance(percentage_loaded, int):
            sin = percentage_loaded * np.reshape(Smax, (len(Smax), 1))
        else:
            sin = np.reshape(percentage_loaded,(len(percentage_loaded),1)) * np.reshape(Smax, (len(Smax), 1))
        sup = np.where(sin>Smax)
        sin[sup[0]] = Smax[sup[0]]
    else:
        sin = np.reshape(Init.sin, (len(Init.sin),1))
    
    sin[0] = BC.edges_bool[0] * BC.edges[0] + (1 - BC.edges_bool[0]) * sin[0]
    sin[-1] = BC.edges_bool[1] * BC.edges[1] + (1 - BC.edges_bool[1]) * sin[-1]
    return sin

def compute_qin(sin, f, k, angle_edges, BC, SD, Init=0):
    if isinstance(Init, int):
        qin = np.dot(Computation.compute_q_from_s(sin), sin)
    
        ## put boundary conditions on Q in the matrix
        qin[0, :] = (1-BC.edges[2])*qin[0, :]
        qin[SD.N_nodes, :] = (1-BC.edges[3])*qin[SD.N_nodes, :]
        qin[0] = BC.edges_bool[2] * BC.edges[2] + [1 - BC.edges_bool[2]] * qin[0]
        qin[-1] = BC.edges_bool[3] * BC.edges[3] + [1 - BC.edges_bool[3]] * qin[-1]
    else:
        qin = np.reshape(Init.Qin, (len(Init.Qin),1))
    return qin
    
def compute_q_sin()