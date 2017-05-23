import numpy as np

"""
    Library defining boundary conditions in xmin and xmax as S, Q or free
    (imposed head,impose flux or free). This properties are stored in two 2
    valued arrays containing boundary types and boundary values. First
    Boundary condition is located downstream (on the river) and the other
    one is located upstream (on the ridge of the hillslope)
    #######################################################################
    INPUTS :
       boundary_type : a two-valued list containing the types corresponding
       to the two boundary_conditions : either S for a constant stock or
       Q for a forced flowrate
       boundary_value : a two-valued list containing the values
       corresponding to each limit. In m/s for Q
    #######################################################################
"""

# Gives an boolean edges matrix from boundary_type
def fixed_edge_matrix_boolean(boundary_type):
    edges_bool = np.zeros(shape=(4, 1))

    if boundary_type[0] == 'S':
        # Start with a fixed charge S for xmin and forced flux Q for xmax
        edges_bool[0] = 1
        edges_bool[2] = 1
    elif boundary_type[0] == 'Q':
        edges_bool[2] = 1

    if boundary_type[1] == 'S':
        edges_bool[1] = 1
        edges_bool[3] = 1
    elif boundary_type[1] == 'Q':
        edges_bool[3] = 1

    return edges_bool

# Gives an edges matrix from boundary_type
def fixed_edge_matrix_values(boundary_type, boundary_value):
    edges = np.zeros(shape=(4, 1))

    if boundary_type[0] == 'S':
        edges[0] = boundary_value[0]
    elif boundary_type[0] == 'Q':
        edges[2] = boundary_value[0]

    if boundary_type[1] == 'S':
        edges[1] = boundary_value[1]
    elif boundary_type[1] == 'Q':
        edges[3] = boundary_value[1]

    return edges
