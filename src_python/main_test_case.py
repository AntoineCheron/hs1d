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
import matplotlib.pyplot as plt
import SaveInitial as SI
import InitialState as IS

###############################################################################
#Flags and input parameters definition
plt.close("all")

#Definition of the flag used to specify the test case
flag = 3
#Definition of folder and name for custom cases (if needed, set flag to 4)
custom_case = "Kerbernez/X_168940_Y_6784091_slopcst"
name_custom = "X_168940_Y_6784091_slopcst"

#Output option for BoussinesqSimulation results
output = 1

#Computation of initial state
initial = 1

#Plot options for vizualisation of results
plot_option = 1
out_put = 0

###############################################################################
#Building output folder simulation results
if not os.path.exists(os.getcwd() + '/simulation_results'):
    os.makedirs(os.getcwd() + '/simulation_results')
    print('Folder for all simulation results saving created')

###############################################################################
#Loading and running model

#Reading the test case and importing data
TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, output=out_put)
print("Test Case files loaded")

#Building output folder for the hillslope
if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
    os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
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
  for i in range(30000):
    rec_init.append(np.mean(TestCase.input.recharge_chronicle))
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
if os.path.exists(os.getcwd() + '/simulation_results/'+ TestCase.name + '/Sin'):
  Init = IS.InitialState(TestCase.name, os.getcwd())
else:
  Init = 0
  perc = 0.5

Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, time_custom=TestCase.input.time, \
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
    out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
    Simulation.output_simu(out_folder)
    print('Simulation results saved')

#Ploting some results
if plot_option !=0:

    plt.figure(num=1)
    plt.plot(Simulation.SR.t[:-1], Simulation.So.recharge_chronicle)
    plt.ylabel('Recharge (m/s)')
    plt.xlabel('t (s)')
    plt.title('Recharge applied to the hillslope as a function of time')
    plt.show()

    plt.figure(num=2)
    plt.plot(Simulation.SD.x_node, Simulation.SD.w_node)
    plt.plot(Simulation.SD.x_edges, Simulation.SD.Hs.w_edges)
    plt.ylabel('Width (m)')
    plt.xlabel('X (m)')
    plt.title('Width of the Hillslope (1) nodes (2) edges')
    plt.show()

    plt.figure(num=3)
    plt.plot(Simulation.SR.t[1:], Simulation.Q_hs)
    plt.ylabel('Flowrate(m3/s)')
    plt.xlabel('t (s)')
    plt.title('Outgoing flowrate as a function of time')
    plt.show()




###############################################################################




