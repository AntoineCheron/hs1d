import numpy as np
from Model import Source as SoU
from Model import BoundaryConditionsUtils as BCU
from Model import SpaceDiscretization as SDU
from Model import InitialConditions as ICU
from Utils import Utils
from Simulation import BoussinesqSimulation as BS
import SimulationResults as SR
from assimulo.problem import Implicit_Problem
from assimulo.solvers.sundials import IDA

def process(Morpho, Geol, Hydro, Init=0, Id='Test'):
    if not isinstance(Morpho.x_custom, int):
        Morpho.nx = np.size(Morpho.x_custom) - 1
        Morpho.xmax = np.max(Morpho.x_custom)
        Morpho.xmin = np.min(Morpho.x_custom)

    #Build Source Class
    So_object = SoU.Source(Hydro.period, Hydro.recharge_type, Hydro.recharge_rate, \
                        Hydro.tmin, Hydro.tmax, Hydro.Nt, Hydro.unit, Hydro.time_custom)
    # Create the source object, containing all the data that we will use
    # accross the program
    SO = Utils.createEmptyObject()
    SO.t = So_object.get_t()
    SO.recharge_chronicle = So_object.get_recharge()
    SO.recharge_type = So_object.get_recharge_type()
    SO.period = So_object.get_period()

    #Build SpaceDiscretization Class
    SD_object = SDU.SpaceDiscretization(Morpho.xmin, Morpho.xmax, Morpho.nx, Morpho.discretization_type, Morpho.x_custom, \
                                     Morpho.angle, Morpho.w, Geol.soil_depth, Geol.k, Geol.f, Morpho.z_custom)

    # Create the SpaceDiscretization object, containing all the data that
    # we will use accross the program
    SD = Utils.createEmptyObject()
    SD.xmin = SD_object.get_xmin()
    SD.xmax = SD_object.get_xmax()
    SD.xcustom = SD_object.get_xcustom()
    SD.discretization_type = SD_object.get_discretization_type()
    SD.angle_node = SD_object.get_angle_node()
    SD.w_node = SD_object.get_w_node()
    SD.soil_depth_node = SD_object.get_soil_depth_node()
    SD.x_node = SD_object.get_x_node()
    SD.dx_node = SD_object.get_dx_node()

    # Retrieve the Hs1D object
    Hs1D = SD.get_hs1d()

    #Set Hydraulic Properties
    k = Hs1D.get_k()
    f = Hs1D.get_f()

    #Build Boundary Conditions
    h0 = f*SD.w_node[0]*SD.soil_depth_node[0]
    Geol.boundary_value[0] = h0

    # Create the SpaceDiscretization object, containing all the data that
    # we will use accross the program
    BC = Utils.createEmptyObject()
    BC.edges_bool = BCU.fixed_edge_matrix_boolean(Geol.boundary_value)
    BC.edges = BCU.fixed_edge_matrix_values(Geol.boundary_type, Geol.boundary_value)

    #Build Initial Conditions
    IC_object = ICU.InitialConditions(Hydro.perc_loaded, f, BC, SD, Hs1D, SO, Init)
    IC = Utils.createEmptyObject()
    IC.sin = IC_object.get_sin()
    IC.qin = IC_object.get_qin()
    IC.q_sin = IC_object.get_q_sin()
    IC.Smax = IC_object.get_Smax()

    #Identification of the hillslope
    Id

    #Mass Matrix of the DAE
    m = BS.compute_mass_matrix()

    # Give all the var we want to Simulate
    return Morpho, Geol, Hydro, Init, Id, SO, SD, k, f, ho, BC, IC, m
