# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:06:10 2017

@author: Quentin Courtois
"""

# The inputs in this file will be the Raw-data from the block
# and the Model
# Geol Morpho Hydro Init Id
import Pre_processing as PreProc
from Simulation import BoussinesqSimulation as BS
from Simulation import SimulationResults as SR
from Utils import Utils
import numpy as np
import sys
import os
pwd = os.getcwd()
sys.path.append(pwd[:-15] + '/src_python/Hillslope1D/process/model_input/')
sys.path.append(pwd[:-15] + '/src_python/Hillslope1D/process/simulation/')
sys.path.append(pwd[:-15] + '/src_python/Hillslope1D/process/tools')
sys.path.append(pwd[:-15] + '/src_python/Hillslope1D/process/tools/boussinesq_tools')
sys.path.append(pwd[:-15] + '/src_python/Hillslope1D/process/tools/file_management')
sys.path.append(pwd[:-15] + '/src_python/Hillslope1D/process/tools/time')
sys.path.append(pwd[:-15] + '/src_python/utils/')

import numpy as np
import LoadCasTest as LCT
import LoadTest as LT
import GeologicInputs as GI
import HydrologicInputs as HI
import MorphologicInputs as MI
from assimulo.problem import Implicit_Problem
from assimulo.solvers.sundials import IDA


#Import Hydro Geol and Morpho
###############################################################################
#Flags and input parameters definition

#Definition of the flag used to specify the test case
flag = 4
#Definition of folder and name for custom cases (if needed, set flag to 4)
custom_case = "straight"
name_custom = "straight"
Id = name_custom

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
#Loading files

#Reading the test case and importing data
TestCase = LCT.LoadCasTest(flag, custom_case, name_custom, output=out_put)
print("Test Case files loaded")

#Building output folder for the hillslope
if not os.path.exists(os.getcwd() + '/../simulation_results/' + TestCase.name):
    os.makedirs(os.getcwd() + '/../simulation_results/' + TestCase.name)
    print('Output folder created for the simulated Hillslope')

z = -1
temp = []
for i in range(len(TestCase.input.time)):
    temp.append(0)
TestCase.input.recharge_chronicle = np.array(temp)
#Building Model geometric inputs
Morpho = MI.MorphologicInputs(nx=len(TestCase.input.x)-1, discretization_type='custom',\
                              x_custom = TestCase.input.x, angle=TestCase.input.i, \
                              z_custom = z, w = TestCase.input.w)

Geol = GI.GeologicInputs(k=TestCase.input.k, f=TestCase.input.f, soil_depth=TestCase.input.soil_depth*5)

perc = 1
Hydro = HI.HydrologicInputs(recharge_rate=TestCase.input.recharge_chronicle, \
                            time_custom=TestCase.input.time, \
                            recharge_type='databased',perc_loaded=perc,unit='sec')

##############################################################################
#Preprocessing to build the Simulation
Init=0
Morpho, Geol, Hydro, Init, Id, SO, SD, Hs1D, ho, BC, IC, m = PreProc.process(Morpho, Geol, Hydro, Init, Id)

# Make SD, Hs1D, SO and BC global to use them into the solver

gO = Utils.createEmptyObject()
gO.SD = SD
gO.Hs1D = Hs1D
gO.SO = SO
gO.BC = BC
gO.m = m

BS.set_global(gO)

###############################################################################
#Resolution of the Simulation
#Building initial state
y0 = np.vstack((IC.sin, IC.qin, IC.q_sin))
y0 = np.squeeze(np.reshape(y0, (len(y0), 1)))

# time range to solve
time_range = SO.t

# Computation of the initial state
yd0 = BS.rhs(time_range[0], y0)

#Build the DAE problem to solve
model = Implicit_Problem(BS.res, y0, yd0, time_range[0])
model.name = 'Boussinesq Simulation, ' + Id
sim = IDA(model)
sim.report_continuously = True # TODO : handle this
sim.verbosity = 10
sim.maxord = 5
sim.maxsteps = 5000

#Solving the DAE using IDA from Assimulo
ncp = len(time_range)
success = 1
try:
    t_res, y_res, yd_res = sim.simulate(time_range[-1], ncp, time_range)
except:
    print("The solver didn't end normally")
    success = 0

print("Success : ", success)
if success == 1:
    SR = SR.SimulationResults(t_res, y_res, SD.N_nodes, Hs1D.N_edges, SD.x_node, Hs1D.x_edges,SO.recharge_chronicle, Hs1D.dx_edges)
