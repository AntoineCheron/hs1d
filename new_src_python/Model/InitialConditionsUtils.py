# -*- coding: utf-8 -*-
"""
Created on Tue May 23 10:27:08 2017

@author: Quentin Courtois
"""
import numpy as np
from Initial import Computation

def compute_Smax(f, SD):

    Smax = f * np.reshape(SD.w_node,(len(SD.w_node),1)) * np.reshape(SD.soil_depth_node,(len(SD.soil_depth_node),1))

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
        sin = np.reshape(Init.Sin, (len(Init.Sin),1))

    sin[0] = BC.edges_bool[0] * BC.edges[0] + (1 - BC.edges_bool[0]) * sin[0]
    sin[-1] = BC.edges_bool[1] * BC.edges[1] + (1 - BC.edges_bool[1]) * sin[-1]

    return sin

def compute_qin(sin, f, k, BC, SD, HS, Init=0):

    if isinstance(Init, int):
        qin = np.dot(-Computation.compute_q_from_s(sin, k, f, SD, HS, BC), sin)
        ## put boundary conditions on Q in the matrix
        qin[0, :] = (1-BC.edges[2])*qin[0, :]
        qin[SD.N_nodes, :] = (1-BC.edges[3])*qin[SD.N_nodes, :]
        qin[0] = BC.edges_bool[2] * BC.edges[2] + [1 - BC.edges_bool[2]] * qin[0]
        qin[-1] = BC.edges_bool[3] * BC.edges[3] + [1 - BC.edges_bool[3]] * qin[-1]
    else:
        qin = np.reshape(Init.Qin, (len(Init.Qin),1))

    return qin

def compute_q_sin(sin, qin, t, f, SD, HS, So, Init=0):

    if isinstance(Init, int):
        q_sin = np.dot(Computation.compute_qs_from_q(np.vstack((sin, \
                        qin, np.zeros((len(qin), 1)))), t[0], f, SD, HS, So), \
                        qin)
    else:
        q_sin = Init.QSin

    temp = np.where(q_sin<0)
    q_sin[temp[0],temp[1]] = 0

    return q_sin
