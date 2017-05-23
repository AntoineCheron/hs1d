# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 09:42:56 2017

@author: Quentin Courtois
"""

import os
import LoadTest as LT
import numpy as np

class BilanCheck(object):
    """
        Class computing the water balance for each simulation in order to check
        the validity of simulations. It is based on computing the difference
        between inputs and outputs of the system

        @param
            name : a list or string containing the name of each simulation
            succes : a list or an integer indicating if the simulation was
                solved successfully (0 or 1)
    """
    def __init__(self, name, success, watershed=""):
        
        #Check if there is only one hillslope
        if watershed == "":
            print('WARNING ! : You didn\'t  refered a watershed\'s name. Creating a generic name base on first hillslope\'s name')
            if isinstance(name, list):
                watershed = 'watershed_' + name[0]
            elif isinstance(name, str):
                watershed = 'watershed_' + name
            else:
                print('Name is not in a valid format')
        #######################################################################
        #Building Simu object which will contains all Simulations
        Simu = []

        #If there are many hillslopes
        if isinstance(name, list):
            for i in range(len(name)):
                if success[i] == 1:
                    folder = os.getcwd() + '/simulation_results/' + watershed + '/' + name[i]
                    Simul_load = LT.LoadTest(folder)
                    Simu.append(Simul_load)
                else:
                    Simu.append(0)
        #If there is one hillslope
        elif isinstance(name, str):
            if success == 1:
                folder = os.getcwd + '/simulation_results/' + name
                Simul_load = LT.LoadtTest(folder)
                Simu.append(Simul_load)
        else:
            print('Wrong type of name')

        #######################################################################
        #Computing the flowrate that goes to river for each simulation
        self.compute_Q_riv(Simu)

        #Computing Stock variation along the hillslope for each time step
        self.compute_storage_variation(Simu)
        #######################################################################
        #Water balance computation
        self.compute_water_balance()
        
        #######################################################################
        #Cumulative water balance computation
        self.compute_cumulative_balance()
        
        #######################################################################
        #Percent_discrepancy computation
        self.compute_percent_discrepancy()
        
        #######################################################################
        #Write water balance computations in a text_file
        self.write_balance(watershed)
        
    def compute_Q_riv(self, Simu):
        """
            Computation of hillslope outflow to the river as a function of time
        """
        #Initializing variables
        Q_riv = np.zeros((len(Simu[0].QS[:,0]),1))
        surf = []

        #Computation for each Simulation
        for i in range(len(Simu)):
            if not isinstance(Simu[i], int):
                temp = Simu[i]
                #Compute the flow that goes in the river as the sum of seepage along the hillslope at each time
                Q_riv = Q_riv + Simu[i].Q_hs
                #Compute the surface of each hillslope
                surf_temp = np.trapz(np.reshape(temp.w_edges, (1, len(temp.w_edges))), \
                                     np.reshape(temp.x_Q, (1, len(temp.x_Q))))
                surf.append(surf_temp)
        
        #Attribute each variable to the class
        self.Q_riv = Q_riv                                 #River flowrate contribution
        self.surface = surf                                #Surface of the hillslopes
        self.t = Simu[0].t_res                             #Times corresponding to teh recharge
        self.recharge = Simu[0].recharge_chronicle         #Recharge chronicle
        
    def compute_water_balance(self):
        """
            Computation of input and output total volumes
        """
        
        #Total volume passed throught the outlet of the river
        Q_riv = np.reshape(self.Q_riv, (1,len(self.Q_riv)))
        t = np.reshape(self.t, (1, len(self.t)))
        self.Q_r = np.trapz(Q_riv, t)
        
        #Total surface of all hillslopes
        self.surface_total = np.sum(self.surface)
        
        #Total volume bringed by the recharge
        self.recharge_vol = np.sum(self.recharge) * self.surface_total
        
    def compute_cumulative_balance(self):
        """
            Comput cumulative water balance for each time
        """
        self.dt = np.diff(np.reshape(self.t,(1,(len(self.t)))))
        cumul_Q = []
        cumul_r = []
        #Initialization
        cumul_Q.append(self.Q_riv[0][0]*self.dt[0][0])
        cumul_r.append(self.recharge[0][0]*self.surface_total*self.dt[0][0])
        for i in range(1,len(self.t)):
            temp_cumul = np.trapz(np.reshape(self.Q_riv[:i+1], (1,len(self.Q_riv[:i+1]))), \
                                    np.reshape(self.t[:i+1], (1, len(self.t[:i+1]))))
            cumul_Q.append(temp_cumul[0])
            temp_cumul = self.surface_total * np.trapz(np.reshape(self.recharge[:i+1], (1,len(self.recharge[:i+1]))), \
                                    np.reshape(self.t[1:i+2], (1, len(self.t[1:i+2]))))
            cumul_r.append(temp_cumul[0])
        
        #Attribute each cumulative
        self.cumul_Q = cumul_Q
        self.cumul_r = cumul_r
    
    def compute_percent_discrepancy(self):
        
        perc_disc = []
        for i in range(len(self.t)):
            perc_disc.append((self.cumul_Q[i] - self.cumul_r[i] + self.dS[i][0])/((self.cumul_Q[i] + self.dS[i][0] + self.cumul_r[i])/2))
        
        self.perc_disc = perc_disc
        
    def compute_storage_variation(self, Simu):
        dS = np.zeros((len(self.t),1))
        ds = []
        for i in range(len(Simu)):
            Sim = Simu[i]
            temp_ds = []
            if not isinstance(Simu[i], int):
                for j in range(len(Sim.t_res)):
                    temp = Sim.S[j,:] - Sim.S[0,:]
                    ds = np.sum(temp)
                    temp_ds.append(ds)
            dS = dS + np.reshape(temp_ds, (len(temp_ds),1))
        self.dS = dS
               
        
    def write_balance(self, watershed):
        """
            Write files containing water balance computation
        """
        #######################################################################
        #Build water_balance and watershed directories
        if not os.path.exists(os.getcwd() + '/water_balance'):
            os.makedirs(os.getcwd() + '/water_balance')
            print('Folder for water balance compuation created')
        if not os.path.exists(os.getcwd() + '/water_balance/' + watershed):
            os.makedirs(os.getcwd() + '/water_balance/' + watershed)
            print('Output folder for the simulated watershed computed')
        
        #######################################################################
        #Save each value in a separated file
        folder = os.getcwd() + '/water_balance/' + watershed + '/'
        
        #Time
        name_file = folder + 't'
        with open(name_file, "wb") as f:
            np.savetxt(f, self.t, fmt='%1.12e', delimiter="\t", newline='\n')
        
        #Q_riv
        name_file = folder + 'Q_riv'
        with open(name_file, "wb") as f:
            np.savetxt(f, self.Q_riv, fmt='%1.12e', delimiter="\t", newline='\n')
        
        #Recharge
        name_file = folder + 'Recharge'
        with open(name_file, "wb") as f:
            np.savetxt(f, self.recharge, fmt='%1.12e', delimiter="\t", newline='\n')        
        
        #Cumulative Q
        name_file = folder + 'Cumulative_Q_vol'
        with open(name_file, "wb") as f:
            np.savetxt(f, np.reshape(self.cumul_Q,(len(self.cumul_Q),1)), fmt='%1.12e', delimiter="\t", newline='\n')
        
        #Cumulative R
        name_file = folder + 'Cumulative_recharge_vol'
        with open(name_file, "wb") as f:
            np.savetxt(f, np.reshape(self.cumul_r,(len(self.cumul_r),1)), fmt='%1.12e', delimiter="\t", newline='\n')

        #Total output volume            
        name_file = folder + 'Output_volume'
        with open(name_file, "wb") as f:
            np.savetxt(f, self.Q_r, fmt='%1.12e', delimiter="\t", newline='\n')        
        
        #Total input volume
        name_file = folder + 'Input_volume'
        with open(name_file, "wb") as f:
            np.savetxt(f, np.array([self.recharge_vol]), fmt='%1.12e', delimiter="\t", newline='\n')

        #Total Surface
        name_file = folder + 'Total_surface'
        with open(name_file, "wb") as f:
            np.savetxt(f, np.array([self.surface_total]), fmt='%1.12e', delimiter="\t", newline='\n')    
            
        #Percent Discrepancy
        name_file = folder + 'percent_discrepancy'
        with open(name_file, "wb") as f:
            np.savetxt(f, self.perc_disc, fmt='%1.12e', delimiter="\t", newline='\n')            
    