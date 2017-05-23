import numpy as np

class Hs1D(object):
    """
        Class containing the hillslope properties. It's an attribute of
        SpaceDiscretization Class
        #######################################################################
        @param
            xmin : minimal coordinate of the hillslope (in meters)
            xmax : maximal coordinate of the hillslope (in meters)
            nx : number of meshes used to describe the hillslope
            angle : slope of the hillslope on each mesh
            w : width of the hillslope on each mesh
            soil_depth : depth of the reservoir on each mesh
            k : hydraulic conductivity on the hillslope (in m/s)
            f : kinematic porosity on the hillslope

        @attributes:
          k : hydraulic conductivity
          f : porosity
          angle_edges : slope on edges
          soil_depth_edges : thickness of the layer on edges
          w_edges : width on edges
        WARNING ! Each property is defined over the edges of the mesh
        #######################################################################
    """

    def __init__(self, nx, angle, w, soil_depth, k, f, z_custom, x):
        #Extracting angle, w, x and soil_depth

        self.x_edges = x

        if isinstance(z_custom,int):
            if isinstance(angle, int) or isinstance(angle, float):
                if angle == 0:
                    self.angle_edges = np.arctan(0.05)*np.ones((nx+1, 1))
                elif angle > 0:
                    self.angle_edges = angle*np.ones((nx+1, 1))
                elif angle < 0:
                    print('Error, angle must be positive. Using absolute value')
                    self.angle_edges = np.abs(angle)*np.ones((nx+1, 1))
            else:
                self.angle_edges = np.reshape(angle, (len(angle), 1))
        else:
            diff_x = np.diff(np.reshape(x,(1,len(x))))
            diff_x = np.append(diff_x,diff_x[-1][-1])
            diff_x = np.reshape(diff_x,(len(diff_x),1))

            diff_z = np.diff(np.reshape(z_custom,(1,len(z_custom))))
            diff_z = np.append(diff_z,diff_z[-1][-1])
            diff_z = np.reshape(diff_z,(len(diff_z),1))
            self.angle_edges = np.arctan(diff_z/diff_x)

        if isinstance(w, int) or isinstance(w, float):
            if w == 0:
                self.w_edges = np.ones((nx+1, 1))*500
            elif w > 0:
                self.w_edges = np.ones((nx+1, 1))*w
            elif w < 0:
                print('Error, width must be positive. Using absolute value')
                self.w_edges = np.ones((nx+1, 1))*np.abs(w)
#            self.w = np.ones((nx+1,1))*100
        else:
            self.w_edges = np.reshape(w, (len(w), 1))

        if isinstance(soil_depth, int) or isinstance(soil_depth, float):
            if soil_depth == 0:
                self.soil_depth_edges = np.ones((1, nx+1))
            elif soil_depth > 0:
                self.soil_depth_edges = np.ones((1, nx+1))*soil_depth
            elif soil_depth < 0:
                print('Error, cells thickness must be a positive value. Using ', \
                      'absolute value')
                self.soil_depth_edges = np.ones((1, nx+1))*np.abs(soil_depth)
        else:
            self.soil_depth_edges = np.reshape(soil_depth, (len(soil_depth),1))

        # Set Hydraulic Parameters
        self.k = k
        self.f = f

    ###########################################################################
    #                               GETTERS                                   #
    ###########################################################################

    def get_w_edges(self):
        return self.w_edges

    def get_soil_depth_edges(self):
        return self.soil_depth_edges

    def get_angle_edges(self):
        return self.angle_edges

    def get_k(self):
        return self.k

    def get_f(self):
        return self.f

    def get_N_edges(self):
        return len(self.x_edges)

    def get_dx_edges(self):
        # length(dx) = Nx
        self.dx_edges = self.x_edges[1:]-self.x_edges[0:-1]
        return self.dx_edges
