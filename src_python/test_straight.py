# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 09:57:09 2017

@author: Quentin Courtois
"""

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
import LoadTest as LT
import SaveInitial as SI
import InitialState as IS
#Flags and input parameters definition
plt.close("all")

#Definition of the flag used to specify the test case
flag = 1

#Definition of folder and name for custom cases (if needed, set flag to 4)
custom_case = ""
name_custom = ""

#Output option for BoussinesqSimulation results
output = 1

#Plot options for vizualisation of results
plot_option = 1
out_put = 1
initial = 1
interp = 0

TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, output=out_put, interp=0)
print("Test Case files loaded")

if not os.path.exists(os.getcwd() + '/simulation_results'):
    os.makedirs(os.getcwd() + '/simulation_results')
    print('Folder for all simulation results saving created')

if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
    os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
    print('Output folder created for the simulated Hillslope')
#Creating initial state
if initial != 0:
    #Building Model inputs
    perc = 0.5
    z=-1
    rec = []
    for j in range(20000):
        rec.append(np.mean(TestCase.input.recharge_chronicle))
    rec = np.reshape(rec,(1, len(rec)))
    toto_toto = rec
    time = np.arange(20000) * 3600
    Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                                  x_custom = TestCase.input.x, angle=TestCase.input.i, \
                                  z_custom = z, w = TestCase.input.w)

    Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

    Hydro = HI.HydrologicInputs(recharge_rate=rec[0], time_custom=time, \
                                recharge_type='databased',perc_loaded=perc, unit= 'hour')

    #Creating the BoussinesqSimulation model using Model Inputs
    Simulation_init = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Id=TestCase.name)
    print("Model built")

    #Running BoussinesqSimulation model
    Simulation_init.implicit_scheme_solver()
#    success = success.append(Simulation.success)
    success = Simulation_init.success
    print("Initial Simulation over")
###############################################################################

###############################################################################
    #Output of the simulation

    #Building output files
    if success == 1:
      save_init = SI.SaveInitial(result=Simulation_init.SR, name=TestCase.name, folder=os.getcwd())
      print("Initial State Saved")
#Reading the test case and importing data
TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, out_put, interp)
print("Test Case files loaded")

z = -1
if os.path.exists(os.getcwd() + '/simulation_results/'+ TestCase.name + '/Sin'):
  Init = IS.InitialState(TestCase.name, os.getcwd())
else:
  Init = 0
  perc = 0.5

#Building Model inputs
Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                              x_custom = TestCase.input.x, angle=TestCase.input.i, \
                              z_custom = z, w = TestCase.input.w)

Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, time_custom=TestCase.input.time, \
                            recharge_type='databased',perc_loaded=perc)

#Creating the BoussinesqSimulation model using Model Inputs
Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Init, Id=TestCase.name)
print("Model built")

#Running BoussinesqSimulation model
Simulation.implicit_scheme_solver()
print("Simulation over\n")
success = Simulation.success

#Output of the simulation

#Building output files
if output != 0:
    if success == 1:
        if not os.path.exists(os.getcwd() + '/simulation_results'):
            os.makedirs(os.getcwd() + '/simulation_results')
            print('Folder for all simulation results saving created')

        if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
            os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
            print('Output folder created for the simulated Hillslope')

        out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
        Simulation.output_simu(out_folder)
        print('Simulation results saved')
    else:
        print('Simulation failed, no results saved')

#Output
#Loading all results
Simu = []
if success == 1:
    folder = os.getcwd() + '/simulation_results/' + TestCase.name
    Simul_load = LT.LoadTest(folder)
    Simu.append(Simul_load)
else:
    Simu.append(0)


#Ploting results
if plot_option != 0:
    Q_riv = 0
    surf = []
    recharge = []
    if success == 1:
        temp = Simu[0]
        Q_riv = Q_riv + np.abs(temp.Q[1:,1]) - np.abs(temp.QS[1:,0]) + \
                np.sum(temp.QS[1:,:],1) + TestCase.input.recharge_chronicle*Simu[0].w_node[0][0]*(Simu[0].x_Q[1][0] - Simu[0].x_Q[0][0])
        surf_temp = np.trapz(np.reshape(temp.w_edges,(1,len(temp.w_edges))), np.reshape(temp.x_Q,(1,len(temp.x_Q))))
        surf.append(surf_temp)
        rec = TestCase.input.recharge_chronicle
    recharge.append(np.trapz(np.reshape(rec,(1,len(rec))), np.reshape(temp.t_res[1:],(1,len(temp.t_res)-1))))

    plt.figure(6)
    plt.plot(Simu[0].t_res[1:], Q_riv, 'b-', Simu[0].t_res, Simulation.Q_hs, 'g--')
    plt.ylabel('Q river (m3/s)')
    plt.xlabel('Time (s)')
    plt.title('Flowrate a the outlet as a function of time')
    plt.grid(True)
    plt.show()



    Q_r = np.trapz(np.reshape(Q_riv,(1,len(Q_riv))), np.reshape(Simu[0].t_res[0:-1], (1, len(Simu[0].t_res[0:-1]))))
    temp_qr = np.trapz(np.reshape(Simulation.Q_hs,(1,len(Simulation.Q_hs))), np.reshape(Simu[0].t_res[0:], (1, len(Simu[0].t_res[0:]))))
    print('Total volume through the outlet : ', float(Q_r), ' m3 over ', int((Simu[0].t_res[-1]+1)/86400), ' days')
    print('Total Surface of the watershed is : ', np.sum(surf), ' m²')
    print('Total Volume In is : ', np.sum(recharge) * np.sum(surf), ' m3 over 35 d')
###############################################################################
