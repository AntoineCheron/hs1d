# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 10:46:35 2017

@author: Quentin Courtois
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:32:18 2017

@author: Quentin Courtois
"""

#Adding all modules path to the main path
import sys
import os
sys.path.append(os.getcwd() + '/Hillslope1D/process/model_input/')
sys.path.append(os.getcwd() + '/Hillslope1D/process/simulation/')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools/boussinesq_tools')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools/file_management')
sys.path.append(os.getcwd() + '/Hillslope1D/process/tools/time')
sys.path.append(os.getcwd() + '/utils/')

import numpy as np
import LoadCasTest as LCT
import GeologicInputs as GI
import HydrologicInputs as HI
import MorphologicInputs as MI
import BoussinesqSimulation as BS
import SaveInitial as SI
import InitialState as IS
import LoadTest as LT

###############################################################################
#Flags and input parameters definition

#Definition of the flag used to specify the test case
flag = 4
#Definition of folder and name for custom cases (if needed, set flag to 4)
custom_case = "straight_slave"
name_custom = "straight_slave"

#Output option for BoussinesqSimulation results
output = 0

#Computation of initial state
initial = 1

out_put = 0

#Recharge
R = 10**-8

#Comparison
comp = 1

###############################################################################
#Building output folder simulation results
if not os.path.exists(os.getcwd() + '/../simulation_results'):
    os.makedirs(os.getcwd() + '/../simulation_results')
    print('Folder for all simulation results saving created')

###############################################################################
#Loading and running model

#Reading the test case and importing data
TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, output=out_put)
print("Test Case files loaded")

#Building output folder for the hillslope
if not os.path.exists(os.getcwd() + '/../simulation_results/' + TestCase.name):
    os.makedirs(os.getcwd() + '/../simulation_results/' + TestCase.name)
    print('Output folder created for the simulated Hillslope')

z = -1

#Building Model geometric inputs
Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x)-1, discretization_type='custom',\
                              x_custom = TestCase.input.x, angle=TestCase.input.i, \
                              z_custom = z, w = TestCase.input.w)

Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

perc = 0.5
if initial != 0:
  #Computing initial recharge applied as the mean value of the applied recharge
  rec_init = []
  for i in range(1000):
    rec_init.append(R)
  t_init = np.arange(len(rec_init))*3600
  Hydro = HI.HydrologicInputs(recharge_rate=rec_init, time_custom=t_init, \
                              recharge_type='databased',perc_loaded=perc)

  #Creating the BoussinesqSimulation model using Model Inputs
  Simulation_init = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Id=TestCase.name)
  print("Initial Model built")

  #Running BoussinesqSimulation model
  Simulation_init.implicit_scheme_solver()
  print("Initial Simulation over")

  #Writing initial State
  if Simulation_init.success != 0:
    save_init = SI.SaveInitial(result=Simulation_init.SR, name=TestCase.name, folder=os.getcwd())
    print("Initial State Saved")
###############################################################################
#Computing real Simulation
if os.path.exists(os.getcwd() + '/../simulation_results/'+ TestCase.name + '/Sin'):
  Init = IS.InitialState(TestCase.name, os.getcwd())
else:
  Init = 0
  perc = 0.5


Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, \
                            time_custom=TestCase.input.time, \
                            recharge_type='databased',perc_loaded=perc)

#Creating the BoussinesqSimulation model using Model Inputs
Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Init, Id=TestCase.name)
print("Model built")

#Running BoussinesqSimulation model
Simulation.implicit_scheme_solver()
print("Simulation over")


###############################################################################
#Output of the simulation

#Building output files
if output != 0:
    out_folder = os.getcwd() + '/../simulation_results/' + TestCase.name
    Simulation.output_simu(out_folder)
    print('Simulation results saved')

###############################################################################
#Comparison with previous simulation
if comp != 0:
    Prev = LT.LoadTest(os.getcwd() +  '/../simulation_results/' + name_custom)

    gap_S = np.sum(Prev.S - Simulation.SR.S) / np.sum((Prev.S + Simulation.SR.S)/2)
    gap_Q = np.sum(Prev.Q - Simulation.SR.Q) / np.sum((Prev.Q + Simulation.SR.Q)/2)
    gap_QS = np.sum(Prev.QS - Simulation.SR.QS) / np.sum((Prev.QS + Simulation.SR.QS)/2)

    print("Error for Stock on all data : " + str(gap_S*100) + " %")
    print("Error for Flowrate on all data : " + str(gap_Q*100) + " %")
    print("Error for Seepage on all data : " + str(gap_QS*100) + " %")

#Test all gaps
try :
    assert gap_S < 0.1
except AssertionError:
    sys.exit("Stock gap is too high")

try :
    assert gap_Q < 0.1
except AssertionError:
    sys.exit("Flow gap is too high")
try :
    assert gap_QS < 0.1
except AssertionError:
    sys.exit("Seepage gap is too high")

###############################################################################
