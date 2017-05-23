# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 10:06:17 2017

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
sys.path.append(os.getcwd() + '/Kerbernez/')

import numpy as np
import LoadCasTest as LCT
import GeologicInputs as GI
import HydrologicInputs as HI
import MorphologicInputs as MI
import BoussinesqSimulation as BS
import InitialState as IS 
import SaveInitial as SI
import BilanCheck as BC

"""Building of a list of hillslopes"""
###############################################################################
#Flags and input parameters definition
watershed_name = 'Kerbernez'
#watershed_name = 'Kerbernez_mono'

#Definition of the flag used to specify the test case
flag = 4
cst = False
out_put = 0
initial_state = 1
#Definition of folder and name for custom cases (if needed, set flag to 4)
src_path = os.getcwd()
if "\\" in src_path:
    src_path = src_path.replace("\\", "/")

watershed_path = src_path[:-10] + 'test_case/matlab/' + watershed_name + '/'
list_dir = os.listdir(watershed_path)

#Check if there is one or more hillslopes
fold_ex = False
for i in range(len(list_dir)):
    if os.path.isdir(watershed_path + '/' + list_dir[i]) is True:
        fold_ex = True

cst_name = []
cst_folder = []
var_name =  []
var_folder = []

#Build a list of hillslopes if there are many hillslopes
if fold_ex is True:
    for i in range(len(list_dir)):
        temp = list_dir[i]
        if temp[-3:] == 'cst':
            cst_name.append(temp)
            cst_folder.append(watershed_name + "/"+temp)
        elif temp[-3:] == 'var':
            var_name.append(temp)
            var_folder.append(watershed_name + "/"+temp)
        else:
            print("Trouble reading watershed's folders")
            
    if cst == True:
        custom_folder = cst_folder
        name_custom = cst_name
        print(len(cst_name),' hillslopes found and loaded')
    else:
        custom_folder = var_folder
        name_custom = var_name
        print(len(var_name),' hillslopes found and loaded')
else:
    name_custom = [watershed_name]
    custom_folder = [watershed_name + '/']
#Output option for BoussinesqSimulation results
output = 1

#Plot options for vizualisation of results
plot_option = 1
6
if not os.path.exists(os.getcwd() + '/simulation_results'):
    os.makedirs(os.getcwd() + '/simulation_results')
    print('Folder for all simulation results saving created')

###############################################################################
"""Compute initial state for all hillslopes"""
#Loading and running model
if initial_state != 0:
    success = np.zeros((len(name_custom),1))
    for i in range(len(name_custom)):
        #Reading the test case and importing data
        TestCase = LCT.LoadCasTest(flag, custom_folder[i], name_custom[i], out_put)
        print("Test Case ", i, " files loaded")
        if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
          os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
          print('Output folder created for the simulated Hillslope')
        perc = 0.2
        z=-1
        rec = []
        for j in range(5000):
            rec.append(np.mean(np.abs(TestCase.input.recharge_chronicle[0:5])))
        rec = np.reshape(rec,(1, len(rec)))
        TestCase.input.recharge_chronicle = rec[0]
        TestCase.input.time = np.arange(5000)
        #Building Model inputs
        Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x), discretization_type='custom',\
                                      x_custom = TestCase.input.x, angle=TestCase.input.i, \
                                      z_custom = z, w = TestCase.input.w)

        Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth)

        Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, time_custom=TestCase.input.time, \
                                    recharge_type='databased',perc_loaded=perc, unit = 'days')

        #Creating the BoussinesqSimulation model using Model Inputs
        Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Id=TestCase.name)
        print("Model built")

        #Running BoussinesqSimulation model
        Simulation.implicit_scheme_solver()
    #    success = success.append(Simulation.success)
        success[i] = Simulation.success
        print("Simulation ", i, " over")

    ###############################################################################

    ###############################################################################
        #Output of the simulation
        if success[i] == 1:
          save_init = SI.SaveInitial(result=Simulation.SR, name=TestCase.name, folder=os.getcwd())
          print("Initial State Saved")
        if i == 0:
          Simu_init = Simulation

        Simulation = 0
print('\n Runing the model\n')
"Loading all hillslopes"
###############################################################################
#Loading and running model
success = np.zeros((len(name_custom),1))
for i in range(len(name_custom)):
    #Reading the test case and importing data
    TestCase = LCT.LoadCasTest(flag, custom_folder[i], name_custom[i], out_put)
    print("Test Case ", i, " files loaded")

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

    Hydro = HI.HydrologicInputs(recharge_rate=np.abs(TestCase.input.recharge_chronicle), time_custom=TestCase.input.time, \
                                recharge_type='databased',perc_loaded=perc, unit='days')

    #Creating the BoussinesqSimulation model using Model Inputs
    Simulation = BS.BoussinesqSimulation(Morpho, Geol, Hydro, Init, Id=TestCase.name)
    print("Model built")

    #Running BoussinesqSimulation model
    Simulation.implicit_scheme_solver()
#    success = success.append(Simulation.success)
    success[i] = Simulation.success
    if i == 0:
        Simu_temp = Simulation
    print("Simulation ", i, " over\n")

###############################################################################

###############################################################################
    #Output of the simulation

    #Building output files
    if output != 0:
        if success[i] == 1:
            if not os.path.exists(os.getcwd() + '/simulation_results'):
                os.makedirs(os.getcwd() + '/simulation_results')
                print('Folder for all simulation results saving created')

            if fold_ex is True:
                if not os.path.exists(os.getcwd() + '/simulation_results/' + watershed_name + '/' + TestCase.name):
                    os.makedirs(os.getcwd() + '/simulation_results/' + watershed_name + '/' + TestCase.name)
                    print('Output folder created for the simulated Hillslope')
                out_folder = os.getcwd() + '/simulation_results/' + watershed_name + '/' + TestCase.name
            else:
                if not os.path.exists(os.getcwd() + '/simulation_results/' + TestCase.name):
                    os.makedirs(os.getcwd() + '/simulation_results/' + TestCase.name)
                    print('Output folder created for the simulated Hillslope')
                out_folder = os.getcwd() + '/simulation_results/' + TestCase.name
                

            
            Simulation.output_simu(out_folder)
            print('Simulation ', i, ' results saved')
        else:
            print('Simulation ', i, ' failed, no results saved')
    Simulation = 0

###############################################################################
#Loading all results and computing Water Balance
Bilan = BC.BilanCheck(name_custom, success, "")
 
###############################################################################
