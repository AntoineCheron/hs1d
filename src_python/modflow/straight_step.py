# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 13:12:19 2017

@author: Quentin Courtois
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 10:29:46 2017

@author: Quentin Courtois
"""
#Add module's path
import sys
import os
sys.path.append(os.getcwd() + '/process/')
sys.path.append(os.getcwd() + '/utils/')

import numpy as np
import ModflowModel2D as MM2D
import WatershedObj as WO
import HydroObj as HO
import WatershedReader as WR
import os
import time
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import copy

#Set Options for the Run
plot_options = 1
output = 1
initial = 1

name = '2017_1_31_9_15_18_X_straight_Y_1_slopcst_Real2';

watershed = WR.WatershedReader(name)
cwd = os.getcwd()

recharge_constant = []
for i in range(4000):
    recharge_constant.append(0.00001)
recharge_true = copy.copy(recharge_constant)
for i in range(4000):
    recharge_true.append(recharge_constant[0]/10)

start_time = time.time()

temporary = np.diff(watershed.coord_x)
temporary = np.unique(temporary)
coord_temporary = np.where(temporary>0)
temporary = temporary[coord_temporary]

hydro = HO.HydroObj(hk=1, sy=0.3, ss=0, soil_depth=2, percentage_loaded=0.1, \
                    recharge_initial = recharge_constant, recharge_true = recharge_true, \
                    time_select=len(recharge_true)-1)

watershed_obj = WO.WatershedObj(name=name, coord_x=watershed.coord_x, coord_y= watershed.coord_y, \
                                elevation=watershed.elevation, bottom=-1, cell_size=temporary[0], \
                                outlet=watershed.outlet)


Mod = MM2D.ModflowModel2D(hydro = hydro, watershed = watershed_obj, output = output, init = initial)
print("------- %s seconds -------" % (time.time() - start_time))



if plot_options != 0:
    Mod.plot_output(name, watershed_obj, hydro)