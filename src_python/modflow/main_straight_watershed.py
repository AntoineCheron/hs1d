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

#Set Options for the Run
plot_options = 1
output = 1
initial = 1

name = '2017_1_31_9_15_18_X_straight_Y_1_slopcst_Real2';

watershed = WR.WatershedReader(name)
cwd = os.getcwd()

rec_temp = np.loadtxt(os.getcwd()+'/Kerbernez/hydrologic.input', skiprows=3)


recharge_true = rec_temp[:,1]*3600
time_rec = rec_temp[:,0]

plt.figure(1)
plt.plot(time_rec, recharge_true)
plt.show()

recharge_constant = np.zeros((1,9000))+np.mean(recharge_true)
recharge_constant = recharge_constant[0]

start_time = time.time()

temporary = np.diff(watershed.coord_x)
temporary = np.unique(temporary)
coord_temporary = np.where(temporary>0)
temporary = temporary[coord_temporary]

hydro = HO.HydroObj(hk=1, sy=0.3, ss=0, soil_depth=2, percentage_loaded=0, \
                    recharge_initial = recharge_constant, recharge_true = recharge_true, \
                    time_select=len(recharge_true)-1)

watershed_obj = WO.WatershedObj(name=name, coord_x=watershed.coord_x, coord_y= watershed.coord_y, \
                                elevation=watershed.elevation, bottom=-1, cell_size=temporary[0], \
                                outlet=watershed.outlet)


Mod = MM2D.ModflowModel2D(hydro = hydro, watershed = watershed_obj, output = output, init = initial)
print("------- %s seconds -------" % (time.time() - start_time))



if plot_options != 0:
    Mod.plot_output(name, watershed_obj, hydro)