# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 11:45:37 2017

@author: Quentin Courtois
"""
import os
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Generate regular watershed
save = 1
plot = 1
name = 'regular_watershed'
X = (np.arange(500))*50
Y = (np.arange(500))*50

coord_o_x = int((len(X)/2)-1)
outlet_x = X[coord_o_x]
outlet_y = 0

Zmax1 = 300
Zmax2 = 50
Zmax3 = 60
Zmin = 0

alpha1 = (Zmax1 - Zmin) / (np.max(Y)-np.min(Y))
beta1 = Zmin
alpha2 = (Zmax1 - Zmax2) / (np.max(Y)-np.min(Y))
beta2 = Zmax2
alpha3 = (Zmax1 - Zmax3) / (np.max(Y)-np.min(Y))
beta3 = Zmax3



Z = np.zeros((len(Y),len(X)))

#Build the 3 main slopes
for j in range(len(Y)):
    Z[j, coord_o_x] = alpha1*Y[j] + beta1
    Z[j, 0] = alpha2*Y[j] + beta2
    Z[j, -1] = alpha3*Y[j] + beta3

for i in range(len(X)):
    Z[-1,:] = Zmax1

alpha_w = np.zeros((len(Y),1))
alpha_e = np.zeros((len(Y),1))
for i in range(len(Y)-1):
    alpha_w[i] = (Z[i, 0] - Z[i, coord_o_x]) / (X[0] - X[coord_o_x])
    alpha_e[i] = (Z[i, -1] - Z[i, coord_o_x]) / (X[-1] - X[coord_o_x])
    for j in range(len(X)):
        if j != coord_o_x and j != 0 and j != len(X)-1:
            if j < coord_o_x:
                Z[i,j] = alpha_w[i]*X[j] + Z[i, 0]
            else:
                Z[i,j] = alpha_e[i]*(X[j]-X[coord_o_x]) + Z[i, coord_o_x]

Elevation = np.zeros((np.size(Z),1))
X_watershed = np.zeros((np.size(Z),1))
Y_watershed = np.zeros((np.size(Z),1))
count = 0

##Set inactives cells
#for i in range(len(Y)):
#    for j in range(len(X)):
#        if i < len(Y)/10 and j < len(X)/10:
#            Z[i,j] = -50
#        elif i < len(Y)/10 and j > 9*len(X)/10:
#            Z[i,j] = -80


for i in range(len(Y)):
    for j in range(len(X)):
        if Z[i,j] > -1:
            X_watershed[count] = X[j]
            Y_watershed[count] = Y[i]
            Elevation[count] = Z[i,j]
            count += 1

if plot > 0:
    fig1 = plt.figure(num=1)
    ax1 = fig1.gca(projection='3d')
    X_new, Y_new = np.meshgrid(X,Y)
    surf = ax1.plot_surface(X_new, Y_new, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig1.colorbar(surf,shrink=0.5, aspect=5)

#save watershed
if save == 1:
   if not os.path.exists(os.getcwd() + "/" + name):
        os.makedirs(os.getcwd() + "/" + name)
   out_folder = os.getcwd() + "/" + name + "/"

   name_file = out_folder + "coord_X_" + name
   with open(name_file, "wb") as f:
       np.savetxt(f, X_watershed, fmt='%1.8e', delimiter="\t", newline='\n')

   name_file = out_folder + "coord_Y_" + name
   with open(name_file, "wb") as f:
       np.savetxt(f, Y_watershed, fmt='%1.8e', delimiter="\t", newline='\n')

   name_file = out_folder + name + "_elevation"
   with open(name_file, "wb") as f:
       np.savetxt(f, Elevation, fmt='%1.8e', delimiter="\t", newline='\n')

   name_file = out_folder + name + "_outlet_coord"
   with open(name_file, "wb") as f:
       np.savetxt(f, [outlet_x, outlet_y], fmt='%1.8e', delimiter="\t", newline='\n')
