# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:27:34 2017

@author: Quentin Courtois
"""

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
import copy
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
plot_option = 0
out_put = 1
initial = 0
interp = 0

TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, output=out_put, interp=0)
print("Test Case files loaded")

rec_init = []
for j in range(4000):
    rec_init.append(0.00001/3600)
#Creating initial state
if initial != 0:
    #Building Model inputs
    perc = 0.5
    z=-1
    rec = copy.copy(rec_init)
    rec = np.reshape(rec,(1, len(rec)))
    time = np.arange(4000) * 3600
    Morpho_init = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                                  x_custom = TestCase.input.x, angle=TestCase.input.i, \
                                  z_custom = z, w = TestCase.input.w)

    Geol_init = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

    Hydro_init = HI.HydrologicInputs(recharge_rate=rec[0], time_custom=time, \
                                recharge_type='databased',perc_loaded=perc, unit= 'hour')

    #Creating the BoussinesqSimulation model using Model Inputs
    Simu_init = BS.BoussinesqSimulation(Morpho_init, Geol_init, Hydro_init, Id=TestCase.name)
    print("Model built")

    #Running BoussinesqSimulation model
    Simu_init.implicit_scheme_solver()
#    success = success.append(Simulation.success)
    success_init = Simu_init.success
    print("Initial Simulation over")
###############################################################################

###############################################################################
    #Output of the simulation

    #Building output files
    if success_init == 1:
        if not os.path.exists(os.getcwd() + '/simulation_results'):
            os.makedirs(os.getcwd() + '/simulation_results')
            print('Folder for all simulation results saving created')

        if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
            os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
            print('Output folder created for the simulated Hillslope')

        out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
        file = out_folder + "/" + TestCase.name + "_init"
        with open(file,"wb") as f:
            np.savetxt(f, Simu_init.SR.S[-1,:], fmt='%1.8e', delimiter="\t", newline='\n')
        file = out_folder + "/" + TestCase.name + "_Smax"
        with open(file,"wb") as f:
            np.savetxt(f, Simu_init.IC.Smax, fmt='%1.8e', delimiter="\t", newline='\n')
        print('Simulation initial state saved\n')
    else:
        print('Simulation failed, no state saved\n')
#run Real Simulation
success = 0
z = -1
if out_put != 0:
    if np.size(TestCase.input.x)-1 == np.size(TestCase.output.storage[:,0]):
        print("using saved initial state")
        perc = TestCase.output.perc_storage[:,0]/100
        z = np.reshape(TestCase.input.z_mod,(len(TestCase.input.z_mod),1))
    else:
        print("using default inital state")
        perc = 0.5
else:
    folder = os.getcwd() + '/simulation_results/' + name_custom + "/" + name_custom
    if os.path.exists(folder + "_init"):
        temp_init = np.loadtxt(folder + "_init")
        temp_Smax = np.loadtxt(folder + "_Smax")
        perc = np.reshape(temp_init, (len(temp_init),1))/np.reshape(temp_Smax, (len(temp_Smax),1))
    else:
        perc = 0.5
#Building Model inputs
z = -1
perc = 0.8
Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                              x_custom = TestCase.input.x, angle=TestCase.input.i, \
                              z_custom = z, w = TestCase.input.w)

Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth, \
                         boundary_type=['A','Q'])

rec_temp = []
#rec_temp =[]
for j in range(4000):
        rec_temp.append(rec_init[0])
for j in range(4000):
        rec_temp.append(rec_init[0]/10)
t_temp = np.arange(12000) * 3600

rec_temp = np.reshape(rec_temp,(1, len(rec_temp)))
rec_temp = rec_temp[0]

Hydro = HI.HydrologicInputs(recharge_rate=rec_temp, time_custom=t_temp, \
                            recharge_type='databased', perc_loaded=perc, unit='hour')

#Creating the BoussinesqSimulation model using Model Inputs
Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Id=TestCase.name)
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
    if success == 1:
        Q_riv = 0
        surf = []
        recharge = []
        if success == 1:
            temp = Simu[0]
            Q_riv = Q_riv + np.abs(temp.Q[:,1]) - np.abs(temp.QS[:,0])*5 + np.sum(temp.QS*5,1)
            surf_temp = np.trapz(np.reshape(temp.w_edges,(1,len(temp.w_edges))), np.reshape(temp.x_Q,(1,len(temp.x_Q))))
            surf.append(surf_temp)
            rec = TestCase.input.recharge_chronicle
        recharge.append(np.trapz(np.reshape(rec,(1,len(rec))), np.reshape(temp.t_res[1:],(1,len(temp.t_res)-1))))

        plt.figure(6)
        plt.plot(Simu[0].t_res, Q_riv)
        plt.ylabel('Q river (m3/s)')
        plt.xlabel('Time (s)')
        plt.title('Flowrate a the outlet as a function of time')
        plt.grid(True)
        plt.show()

        Q_r = np.trapz(np.reshape(Q_riv,(1,len(Q_riv))), np.reshape(Simu[0].t_res, (1, len(Simu[0].t_res))))
        print('Total volume through the outlet : ', float(Q_r), ' m3 over ', int((Simu[0].t_res[-1]+1)/86400), ' days')
        print('Total Surface of the watershed is : ', np.sum(surf), ' m²')
        print('Total Volume In is : ', np.sum(recharge) * np.sum(surf), ' m3 over 35 d')
    ###############################################################################
