# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:02:09 2017

@author: Quentin Courtois
"""

import numpy as np


class InitialState(object):
  """
    Class loading precomputed initial state S, Q, QS
    @param
      name : name of the hillslope computed
      folder : folder containing the folder simulation_results

    @attributes
      Sin : initial stock along the hillslope (on nodes)
      Qin : initial flowrate along the hillslope (on edges)
      QSin : initial seepage along the hillslope (on nodes)
  """

  def __init__(self, name, folder):
    file_S = folder + '/simulation_results/' + name + '/Sin'
    file_Q = folder + '/simulation_results/' + name + '/Qin'
    file_QS = folder + '/simulation_results/' + name + '/QSin'

    self.Sin = self.load_initial_state(file_S)
    self.Qin = self.load_initial_state(file_Q)
    self.QSin = self.load_initial_state(file_QS)

  def load_initial_state(self,file):
    temporary = np.loadtxt(file)
    out = temporary[-1,:]
    return out