# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 10:05:30 2017

@author: Quentin Courtois
"""

import numpy as np
from scipy.interpolate import UnivariateSpline

class GeometricTool(object):
    
    def __init__(self, x, w, angle, soil_depth, k=4, s=0.9):
        #Initialization
        self.x_edges = x
        self.w_edges = w
        self.soil_depth_edges = soil_depth
        self.angle_edges = angle
        self.k = k
        self.s = s
        
        #Compute x_node
        self.x_node = self.compute_x_node()
    
        
    def resample_edges(self):   
        
        #Compute w_node
        smooth_width_function = UnivariateSpline(self.x_edges,self.w_edges, k=self.k, s=self.s)
        self.w_node = smooth_width_function(self.x_node)
        
        #Compute angle_node
        smooth_slope_function = UnivariateSpline(self.x_edges, self.angle_edges, s=0.9)
        self.angle_node = smooth_slope_function(self.x_node)
        self.soil_depth_node = np.interp(np.squeeze(self.x_node),np.squeeze(self.x_edges), self.soil_depth_edges[0, :])
        
    def compute_x_node(self):
        x_node = (self.x_edges[1:] + self.x_edges[0:-1])/2
        return x_node