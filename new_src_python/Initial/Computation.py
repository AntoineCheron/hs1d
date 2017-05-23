# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:16:46 2017

@author: Quentin Courtois
"""

import numpy as np

def compute_q_from_s(y, k, f, SD, HS, BC):
    """
        Compute the flow rate in each cell (over the edges) using conversion
        matrices and stock value
    """
    tempo = y[:SD.N_nodes]/(f*np.reshape(SD.w_node, (SD.N_nodes, 1)))
    q_from_s = (k/f)*(np.cos(HS.angle_edges)*np.dot(SD.b, tempo) \
                + np.sin(HS.angle_edges))
    q_from_s = np.diag(q_from_s[:, 0])
    q_from_s = np.dot(q_from_s, SD.omega)

    ## put boundary conditions on Q in the matrix
    q_from_s[0, :] = (1-BC.edges_bool[2])*q_from_s[0, :]
    q_from_s[SD.N_nodes, :] = (1-BC.edges_bool[3])*q_from_s[SD.N_nodes, :]
    
    return q_from_s

def compute_qs_from_q(y, t, alpha, SD):
    """
        Compute Seepage in each cell (over nodes) based  on flow rate and
        conversion matrices
    """
    alpha_complementar = np.diag(1 - alpha[:, 0])
    qs_from_q = np.dot(-alpha_complementar, SD.a)

    return qs_from_q

def compute_alpha(y, t, f, SD, HS, recharge):
    """
        compute a matrix defining variations over time. Used to compute Q,
        S and QS
    """
    SD.w_node = np.reshape(SD.w_node, (len(SD.w_node), 1))
    SD.soil_depth_node = np.reshape(SD.soil_depth_node, \
                                         (len(SD.soil_depth_node), 1))
    xt = y[0:SD.N_nodes]/(f* SD.soil_depth_node * SD.w_node)
    z = thresholdfunction(xt)
    test_deriv = test_derivative(y, t, SD, HS, recharge)
    alpha =  z * test_deriv + (1 - test_deriv)
    return alpha

def compute_dsdt_from_q(y, t, alpha, beta, SD):
    """
        Compute stock variation between two time steps based on flow rate and
        conversion matrices
    """
    dsdt_from_q = np.dot(np.dot(-np.diag(beta), np.diag(alpha)),SD.a)

    return dsdt_from_q


def test_derivative(y, t, SD, HS, recharge):
    """
        Test to determine if stock is still positive and seepage is occuring
        or not
    """
    
    recharge_rate_spatialized = recharge*np.reshape(SD.w_node, (len(SD.w_node), 1))

    test_deriv = np.dot(-SD.a, y[SD.N_nodes: SD.N_nodes + HS.N_edges]) + recharge_rate_spatialized
    test_deriv = test_deriv >= 0

    return test_deriv

def compute_recharge_rate(t, t_chronicle, recharge_type, period, recharge_chronicle):
    if period == 'inf':
        recharge = recharge_chronicle

    elif period == 0:
        recharge = np.interp(t, t_chronicle, recharge_chronicle)

    elif period == -1:
        recharge = np.interp(t, t_chronicle, recharge_chronicle)

    elif recharge_type == 'periodical':
        recharge = recharge_chronicle*(1+np.cos(2*np.pi*(t/period)))
    elif recharge_type == 'square':
        Int_ = np.floor(t/period)
        Odd_rest = Int_%2
        recharge = recharge_chronicle*Odd_rest
        rem_stiff = t%(2*period)
        bool1 = rem_stiff < 60
        if bool1 is True:
            recharge = recharge_chronicle*(1-(rem_stiff)/60)
        rem_stiff = (t + period)%(2*period)
        bool1 = rem_stiff < 60
        if bool1 is True:
            recharge = recharge_chronicle*((rem_stiff)/60)
    return recharge