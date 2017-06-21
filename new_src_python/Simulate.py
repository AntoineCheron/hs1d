# The inputs in this file will be the Raw-data from the block
# and the Model
# Geol Morpho Hydro Init Id
import Pre_processing as PreProc
from Simulation import BoussinesqSimulation as BS
from Simulation import SimulationResults as SR
from Utils import Utils
import numpy as np



#Preprocessing to build the Simulation
Morpho, Geol, Hydro, Init, Id, SO, SD, Hs1D, ho, BC, IC, m = PreProc.process(Morpho, Geol, Hydro, Init, Id)

# Make SD, Hs1D, SO and BC global to use them into the solver
globalObject = Utils.createEmptyObject()
globalObject.SD = SD
globalObject.Hs1D = Hs1D
globalObject.SO = SO
globalObject.BC = BC
global globalObject

###############################################################################
#Resolution of the Simulation
#Building initial state
y0 = np.vstack((IC.sin, IC.qin, IC.q_sin))
y0 = np.squeeze(np.reshape(y0, (len(y0), 1)))

# time range to solve
time_range = SO.t

# Computation of the initial state
yd0 = BS.res(time_range[0], y0)

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
    SR = SR.SimulationResults(t_res, y_res, SD.N_nodes, Hs1D.N_edges, SD.x_node, Hs1D.x_edges,So.recharge_chronicle, Hs1D.dx_edges)
