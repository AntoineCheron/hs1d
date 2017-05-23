# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:18:56 2017

@author: Quentin Courtois
"""
import numpy as np

class SaveInitial(object):

  """
    Class writing output of the simulation for initial state computations
    @param
      name : name of the hillslope computed
      folder : folder containing the folder simulation_results

    @attributes
      herited from BoussinesqSimulation.SimulationResults
  """


  def __init__(self, result, name, folder):
    fold = folder + '/simulation_results/' + name + '/'

    self.write_initial_state(fold, 'Sin', result.S)
    self.write_initial_state(fold, 'Qin', result.Q)
    self.write_initial_state(fold, 'QSin', result.QS)


  def write_initial_state(self, folder, var_name, var):
    file = folder + var_name
    with open(file,"wb") as f:
      np.savetxt(f, var, fmt='%1.12e', delimiter="\t", newline='\n')