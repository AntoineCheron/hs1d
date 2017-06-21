import numpy as np
from scipy.interpolate import UnivariateSpline
from Model import Hs1D as Hs
import copy

class SpaceDiscretization(object):
    """
        Class defining the discretization in space of the hillslope (linear, log, etc....)
        It's detirmined using X minimal and maximal values, soil depth, angle and width
        of the hillslope for each mesh

        @param :
          - xmin : minimal coordinate of the hillslope
          - xmax : maximal coordinate of the hillslope
          - nx : number of cells to define (if x_custom == -1)
          - discretization_type : type of discretization used : 'linear', 'logarithmic',
              'square' or 'custom'
          - x_custom : -1 or vector of coordinates if discretization_type == 'custom'
          - angle : float or vector of values corresponding to the slope
          - w : float or vector of values corresponding to the width
          - soil_depth : float or vector of values corresponding to the thickness of the layer
          - k : hydraulic conductivity
          - f: porosity
          - z_custom : -1 or a vector altitude of each cell (used to compute the slope
               if x_custom is a vector)

        @attributes:
          - xmin : minimal coordinate of the hillslope
          - xmax : maximal coordinate of the hillslope
          - discretization : type of discretization used
          - xcustom : -1 or vector of coordinates if discretization_type == 'custom'
          - b : computation matrix
          - a : computation matrix
          - omega : computation matrix
          - omega2 : computation matrix
          - x_node : coordinates of cells nodes
          - angle_node : slope of cells on nodes
          - w_node : width of cells on nodes
          - N_node : number of cells
          - soil_depth_node : layer thickness on nodes
          - dx_node : distance between two successive nodes
          - x_edges : coordinates of cells edges
          - N_edges : number of cells edges
          - dx_edges : distance between two successive edges
          - Hs : class Hs1D
    """

    def __init__(self, xmin=0, xmax=100, nx=100, discretization_type='linear', x_custom=-1, \
                 angle=0, w=0, soil_depth=0, k=1/3600, f=0.3, z_custom=-1):
        self.xmin = xmin
        self.xmax = xmax
        self.xcustom = x_custom
        self.discretization_type = discretization_type
        self.x_edges, self.N_edges, self.xmin, self.xmax = self.space_discretization(discretization_type, xmin, xmax, x_custom)
        self.Hs = Hs.Hs1D(nx, angle, w, soil_depth, k, f, z_custom, self.x_edges)
        self.x_node = self.get_x_node()
        self.a, self.b, self.omega, self.omega2 = self.get_matrix_properties()

    ###########################################################################
    #                               GETTERS                                   #
    ###########################################################################

    def get_xmin(self):
        return self.xmin

    def get_xmax(self):
        return self.xmax

    def get_xcustom(self):
        return self.xcustom

    def get_discretization_type(self):
        return self.discretization_type

    def get_hs1d(self):
        return self.Hs

    def get_a(self):
        return self.a

    def get_b(self):
        return self.b

    def get_omega(self):
        return self.omega

    def get_omega2(self):
        return self.omega2

    def get_angle_node(self):
        w_node, angle_node, soil_depth_node = self.resample_hs1D_spatial_variables()
        return angle_node

    def get_w_node(self):
        w_node, angle_node, soil_depth_node = self.resample_hs1D_spatial_variables()
        return w_node

    def get_soil_depth_node(self):
        w_node, angle_node, soil_depth_node = self.resample_hs1D_spatial_variables()
        return soil_depth_node

    def get_x_node(self):
        self.x_node = (self.x_edges[1:] + self.x_edges[0:-1])/2
        return self.x_node

    def get_dx_node(self):
        # length(dxS) = Nx - 1
        self.dx_node = self.x_node[1:]-self.x_node[0:-1]
        return self.dx_node
    
    def get_dx_edges(self):
        self.dx_edges = self.x_edges[1:] - self.x_edges[0:-1]
        return self.dx_edges
    
    def get_N_nodes(self):
        return len(self.x_node)
    
    def get_N_edges(self):
        return len(self.x_edges)

    ###########################################################################
    #                           CUSTOM METHODS                                #
    ###########################################################################


    def space_discretization(self, discretization, xmin, xmax, xcustom):
        # Create the variable that will be returned at the end
        x_edges = None
        N_nodes = None
        N_edges = None

        if discretization == 'linear':
            x_edges = np.arange(xmin, xmax, (xmax - xmin) / N_edges)
        elif discretization == 'logarithmic':
            if xmin == 0:
                log_xmin = -2
                xmin = 10**(log_xmin)
            else:
                log_xmin = np.log10(xmin)
            log_xmax = np.log10(xmax)
            x_edges = np.logspace(log_xmin, log_xmax, N_edges)
        elif discretization == 'custom':
            x_edges = xcustom
            N_edges = len(x_edges)
            N_nodes = N_edges - 1
            xmin = min(x_edges)
            xmax = max(x_edges)
        elif discretization == 'square':
            x_edges = np.arange(xmin**0.5,xmax**0.5, \
                                     (xmax**0.5 - xmin**0.5)/N_edges)
            x_edges = x_edges**2
        else:
            print('no discretization type corresponding to ', discretization, '\n')

        return x_edges, N_edges, xmin, xmax

    def resample_hs1D_spatial_variables(self):
        soil_depth = self.Hs.get_soil_depth_edges()
        w = self.Hs.get_w_edges()
        angle = self.Hs.get_angle_edges()
        x_node = self.get_x_node()
         #UnivariateSpline does the same work as fit ('smoothingspline') in matlab
        smooth_width_function = UnivariateSpline(self.x_edges, w, k=4, s=0.9)
        w_node = smooth_width_function(x_node)
        smooth_slope_function = UnivariateSpline(self.x_edges, angle, s=0.9)
        angle_node = smooth_slope_function(x_node)
        soil_depth_node = np.interp(np.squeeze(x_node),np.squeeze(self.x_edges), soil_depth[0, :])
        return w_node, angle_node, soil_depth_node

    def get_matrix_properties(self):
        a = self.first_derivative_downstream()
        b = self.first_derivative_upstream()
        omega = self.weight_matrix(b)
        omega2 = self.weight_matrix_bis(b)
        return (a, b, omega, omega2)



        ##########################################################################
        ###### stencil matrix(for basic derivation centered, downstream, upstream)
        ###### + mean matrix (Omega)######
        ##########################################################################

    def first_derivative_upstream(self):
        #######################################
        #         [ 0   ...   0   0   ...   0 ]
        #         [-1    1    0   0   ...   0 ]
        #    B =  [ 0  - 1    1   0   ...   0 ]
        #         [... ...   ... ...  ...  ...]
        #         [ 0  ...   ...  0   - 1   1 ]
        #         [ 0  ...    0   0   ...   0 ]
        #size B (Nx + 1) x(Nx)
        #######################################

        a = self.first_derivative_downstream()
        a = copy.copy(a[:-1, :-1])
        a[a > 0] = 1
        a[a < 0] = -1
        temp = np.zeros((1, (len(self.x_edges)-1)))
        b = np.append(temp, a, axis=0) #+ np.zeros(((len(self.x) - 1),1))
        b = np.append(b, temp, axis=0)
#        #B = [zeros(1, length(obj.x) - 1); A; zeros(1, length(obj.x) - 1)];
        dx_node = self.get_dx_node(),
        dx_node = np.hstack((0, np.squeeze(dx_node), 0))
        b = np.dot(np.diag(1/dx_node), b)
        b[np.isnan(b)] = 0
        # self.B = np.eye(1 / self.dxS) * self.B
        #self.b = csr_matrix(self.b)
        return b

    def first_derivative_downstream(self):
        #########################################
        #       [-1   1   0    0    0   ...   0 ]
        # A =   [ 0 - 1   1    0    0   ...   0 ]
        #       [... ... ...  ...  ...  ...  ...]
        #       [ 0  ... ...  ...   0  - 1    1 ]
        #       size A = (Nx) x(Nx + 1)
        #########################################
        dx_edges = self.get_dx_edges()
        a = -np.eye(len(self.x_edges)-1)
        temp = np.zeros(((len(self.x_edges)-1), 1))
        a = np.hstack((a, temp))
        a1 = np.eye(len(self.x_edges) - 1)
        a2 = np.zeros((len(self.x_edges) - 1, 1))
        abis = np.hstack((a2, a1))
        a = a + abis
        a = np.dot(np.diag(np.squeeze(1/self.dx_edges)), a)
#        print(np.size(self.a,0),'  ',np.size(self.a,1))
        # self.A = (np.eye(np.ones(((len(self.x)-1,1)))))*self.A
        #self.A = csr_matrix(self.A)

        return a

    def first_derivative_centered(self):
        ######################################################
        #           [ 0      1       0    0   ...           0]
        # C =       [ -1     0       1    0   ...           0]
        #           [ ...   ...     ...                     1]
        #           [ 0     ...     ...        0    - 1     0]
        ######################################################

        c = -np.diag((np.ones((len(self.x_edges) - 1, 1)), -1))
        c = c + np.diag(np.ones((len(self.x_edges) - 1, 1)), 1)
        c = np.diag(1 / self.dx_edges) * c

        return c

    def weight_matrix(self,b):
        omega = copy.copy(b)
        omega[(omega > 0)] = 0.5
        omega[(omega < 0)] = 0.5
        return omega

    def weight_matrix_bis(self,b):
        omega2 = copy.copy(b)
        omega2[omega2 < 0] = 0
        return omega2
