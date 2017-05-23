import numpy as np
from Model import Source as So
from Model import BoundaryConditionsUtils as BC
from Model import SpaceDiscretization as SD
from Model import InitialConditions as IC
from Initial import Computation
import SimulationResults as SR
from assimulo.problem import Implicit_Problem
from assimulo.solvers.sundials import IDA

"""
###########################################################################
Library called by BoussinesqSimulation

- discretization : space_discretization object contains also the spatial
  properties of your hillslopes (x, width function w and soil_depth)
      named : SD
- source_terms : source object containing the time properties and the
  recharge related to it
      named : So
- boundary_conditions : boundary object containing the time of boundary on
  Morpho.xmin & Morpho.xmax (Q fixed, S fixed or free condtions)
      named : BC
- initial_conditions : initial conditions for S, Q & QS
      named : IC
- Integration results are contained in self.t_res, self.y_res, self.yd_res
      named : SR

###########################################################################
INPUT format and units:

  @param:
    Inputs contained in 3 different classes
    - Geol : k, f and soil_depth (parameters)
    - Morpho : nx, xmin, xmax, discretization_type, x_custom, angle,
               z_custom, w (geometry : spatial aspect)
    - Hydro : tmin, tmax, Nt, unit, recharge_rate, time_custom, period,
              recharge_type, perc_loaded (temporal aspect)
	- Init : Sin, Qin, Qsin if a previous initial state was computed

        WARNING ! ALL CALCULATIONS ARE DONE IN M/S

- k : hydraulic conductiviy (m/s) (ranges from e-8 to e-2 m/s)
- f : porosity (in %) (ranges from e-5 to 0.3)
- soil_depth : thickness of the layer (m) (ranges from 2 to 50)
- period : used only if recharge_type is not data based : in the unit as t
- recharge_type : either 'square','periodical','steady','random','databased'
- recharge_rate : only if not data_based : in mm/day

- tmin : minimal time of the serie
- tmax : maximal time of the serie
- Nt : Number of time steps in the time serie
- unit : time unit of the serie

- Geol.boundary_type : type of the two boundary conditions, first one corresponds
  to the boundary at the level of the river, second at the top of the
  hillslope
      Example : ['S','Q' ] with S : imposed level and Q : imposed flow rate
- boundary_values : values corresponding to the two boundary conditions
      Example : [0,0] in m and in m/sec

- Morpho.xmin : minimal coordinate of the hillslope (m)
- Morpho.xmax : maximal coordinate of the hillslope (m)
- discretization space along the hillslope (int) : default is 100
- Morpho.discretization_type : type of the segmentation of the hillslope :
    'linear', 'logarithmic','square','custom'
- Morpho.x_custom : used only for custom segmentation : must contain edges coordinates
- Morpho.angle : Morpho.angle of the hillslope (rad)
- w : width of the hillslope for each x
- soil_depth : thickness of the media for each x
- k : the kincematic porosity (in % : 30% = 0.3)
- f : the hydraulic conductivity (m/second)

- Id : the identification of the hillslope
- percentage_loaded : initial charge in the hillslope (a percentage of the
  total capacity as 100% = 1)
###########################################################################
"""

def compute_c(self, y, t, SD, k, f, Hs1D, SO, BC):
    """
        Compute a part of the DAE base on S, Q and QS values
    """
    ############################################
    #       [ 0      - alpha * beta * A      0 ]
    # C =   [P(y)       I                    0 ]
    #       [ 0      - (1 - alpha) * A     - I ]
    ############################################

    c = np.zeros((2 * SD.N_nodes + SD.N_edges, 2 * SD.N_nodes + SD.N_edges))

    dsdt_from_q = Computation.compute_dsdt_from_q(y, t, f, SD, Hs1D, SO.recharge_chronicle)
    q_from_s = Computation.compute_q_from_s(y, k, f, SD, Hs1D, BC)
    qs_from_q = Computation.compute_qs_from_q(y, t, f, SD, Hs1D, SO.recharge_chronicle)

        # first row block
    c[0:SD.N_nodes, 0:SD.N_nodes]                                          = 0
    c[0:SD.N_nodes, SD.N_nodes:(SD.N_nodes + SD.N_edges)]                  = dsdt_from_q
    c[0:SD.N_nodes, (SD.N_nodes + SD.N_edges):(2*SD.N_nodes + SD.N_edges)] = 0

        # second row block
    c[SD.N_nodes :(SD.N_nodes + SD.N_edges), 0:SD.N_nodes]                                        = q_from_s
    c[SD.N_nodes :(SD.N_nodes + SD.N_edges), SD.N_nodes:SD.N_nodes + SD.N_edges]                  = np.eye(SD.N_edges)
    c[SD.N_nodes :(SD.N_nodes + SD.N_edges), (SD.N_nodes + SD.N_edges):2*SD.N_nodes + SD.N_edges] = 0

        # third row block
    c[(SD.N_nodes + SD.N_edges):2*SD.N_nodes + SD.N_edges, 0:SD.N_nodes]                                        = 0
    c[(SD.N_nodes + SD.N_edges):2*SD.N_nodes + SD.N_edges, (SD.N_nodes):SD.N_nodes + SD.N_edges]                = qs_from_q
    c[(SD.N_nodes + SD.N_edges):2*SD.N_nodes + SD.N_edges, (SD.N_nodes + SD.N_edges):2*SD.N_nodes + SD.N_edges] = -np.eye(SD.N_nodes)

    # Introduce boundary conditions in the matrix
    c[0, :] = (1-int(BC.edges_bool[0]))*c[0, :]
    c[SD.N_nodes-1, :] = (1-BC.edges_bool[1])*c[SD.N_nodes-1, :]

    return c


def compute_source_terms(self, y, t):
    """
        compute the recharge for a time step on each cell based on recharge
        defined by user
    """
    #Matrix describing source term (recharge)
    N_edges = self.SD.N_edges
    N_nodes = self.SD.N_nodes
    d = np.zeros((2*N_nodes + N_edges, 1))
    alpha = self.compute_alpha(y, t)
    beta = self.compute_beta(y,t)
    Recharge_rate = self.So.compute_recharge_rate(t)
    Recharge_rate_spatialized = Recharge_rate * self.SD.w_node

    d[0:N_nodes] = beta * alpha * Recharge_rate_spatialized
    d[(N_nodes + N_edges):2*N_nodes + N_edges] = (1 - alpha) * Recharge_rate_spatialized


    edges_bool = self.BC.fixed_edge_matrix_boolean()
    edges = self.BC.fixed_edge_matrix_values()

    d[0, :] = (1-edges_bool[0])*d[0, :]
    d[N_nodes-1, :] = (1-edges_bool[1])*d[N_nodes-1, :]
    d[N_nodes, :] = (1-edges_bool[2])*d[N_nodes, :]
    d[N_nodes + N_edges-1, :] = (1-edges_bool[3])*d[N_nodes + N_edges - 1, :]
#       Set D(edges) at the fixes value (for Q)
    d[N_nodes, :] = -edges_bool[2] * edges[2]
    d[N_nodes + N_edges-1, :] = -edges_bool[3] * edges[3]
    d[d == -0] = 0

    return d


def rhs(t, y):
    """
        Compute differential equation dy/dt = C*dy + d
    """
    ##########################################
    #Differential equation : dy/dt = c*dy + d#
    ##########################################
    y = np.reshape(y, (len(y), 1))
    c = bouss_obj.compute_c(y, t)
    d = bouss_obj.compute_source_terms(y, t)
    dy = np.dot(c, y)+d
    dy = np.squeeze(np.reshape(dy, (np.size(dy), 1)))
    return dy

def res(t, y, yd):
    """
        Compute algebraic differential equation m * dy/dt = C*dy + d
    """
    #####################################################
    #function describing the DAE as MM * dy/dt = c*y + d#
    #####################################################
    yd = np.reshape(yd, (len(yd), 1))
    y = np.reshape(y, (len(y), 1))
    dy_new = np.dot(bouss_obj.compute_c(y, t), y) + bouss_obj.compute_source_terms(y, t)
    res = np.dot(bouss_obj.m, yd) - np.reshape(dy_new, (len(dy_new), 1))
    res = np.squeeze(np.reshape(res, (len(res), 1)))
    return res

def implicit_scheme_solver(self):
    """
        Resolution of the DAE using implicit problem solver DAE
    """
    #Building initial state
    y0 = np.vstack((self.IC.sin, self.IC.qin, self.IC.q_sin))
    y0 = np.squeeze(np.reshape(y0, (len(y0), 1)))

    f = BoussinesqSimulation.rhs
    # time range to solve
    t = self.So.TP.t

    global bouss_obj
    bouss_obj = self
    # Computation of the initial state
    yd0 = f(t[0], y0)

    #Build the DAE problem to solve
    self.model = Implicit_Problem(BoussinesqSimulation.res, y0, yd0, t[0])
    self.model.name = 'Boussinesq Simulation'
    self.sim = IDA(self.model)
    self.sim.report_continuously = True
    self.sim.verbosity = 10
    self.sim.maxord = 5
    self.sim.maxsteps = 5000

    #Solving the DAE using IDA from Assimulo
    ncp = len(t)
    ncp_list = t
    self.success = 1
    try:
        t_res, y_res, yd_res = self.sim.simulate(t[-1], ncp, ncp_list)
    except:
        print("The solver didn't end normally")
        self.success = 0

    print("Success : ", self.success)
    bouss_obj = 0
    if self.success == 1:
        self.SR = SR.SimulationResults(t_res, y_res, self.SD.N_nodes, \
                                   self.SD.N_edges, self.SD.x_node, self.SD.x_edges,\
                                   self.So.recharge_chronicle, self.SD.dx_edges)

def compute_mass_matrix(self):
    """
        Compute the matrix m of the DAE
    """
    #Mass matrix used as a multiplier for dy/dt in DAE
    N_edges = self.SD.N_edges
    N_nodes = self.SD.N_nodes
    n_unknowns = 2*N_nodes + N_edges
    m = np.zeros((n_unknowns, n_unknowns))
    m[0:N_nodes, 0:N_nodes] = np.eye(N_nodes)
    return m

def output_simu(self, folder):
    """
        write Simulations Results (Q,S,QS and x_Q,x_S,t_res) in  .txt files
        (delimiter : tab) in the current working directory
    """
    x_S = self.SR.x_node
    t = self.SR.t
    x_Q = self.SD.x_edges

    S_sol = self.SR.S
    Q_sol = self.SR.Q
    QS_sol = self.SR.QS
    #Computation of the outgoing flow of the hillslope (m3/s)
    self.Q_hs = np.sum(self.SR.QS,1)

    #Save Stock integration results
    name_file = folder + "/S"
    with open(name_file, "wb") as f:
        np.savetxt(f, S_sol, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save flux integration results
    name_file = folder + "/Q"
    with open(name_file, "wb") as f:
        np.savetxt(f, Q_sol, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save Seepage integration results
    name_file = folder + "/QS"
    with open(name_file, "wb") as f:
        np.savetxt(f, QS_sol, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save time requested points during integration
    name_file = folder +"/t_res"
    with open(name_file, "wb") as f:
        np.savetxt(f, t, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save mesh's centers coordinates (Morpho.nx)
    name_file = folder +"/x_S"
    with open(name_file, "wb") as f:
        np.savetxt(f, x_S, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save mesh's edges coordinates (Morpho.nx+1)
    name_file = folder +"/x_Q"
    with open(name_file, "wb") as f:
        np.savetxt(f, x_Q, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save recharge chronicle
    name_file = folder +"/recharge"
    with open(name_file, "wb") as f:
        np.savetxt(f, self.So.recharge_chronicle, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save recharge's time valeus
    name_file = folder +"/time"
    with open(name_file, "wb") as f:
        np.savetxt(f, self.So.TP.t, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save maximal storga of each cell
    name_file = folder +"/Smax"
    with open(name_file, "wb") as f:
        np.savetxt(f, self.IC.Smax, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save width of hillslope's cells
    name_file = folder +"/width_S"
    with open(name_file, "wb") as f:
        np.savetxt(f, self.SD.w_node, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save wodth of hillslope's cells edges
    name_file = folder +"/width_Q"
    with open(name_file, "wb") as f:
        np.savetxt(f, self.SD.Hs.w_edges, fmt='%1.12e', delimiter="\t", newline='\n')
    #Save flow of the hillslope
    name_file = folder +"/Q_hillslope"
    with open(name_file, "wb") as f:
        np.savetxt(f, self.Q_hs, fmt='%1.12e', delimiter="\t", newline='\n')
