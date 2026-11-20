import numpy as np
#cell ids
VACUUM = 0
WATER = 1
FUEL = 2
#other constants
ASSEMBLY_PITCH = 21

#square mesh
def make_square_mesh(n_assemblies):
    assembly_map = np.full((n_assemblies, n_assemblies), FUEL)
    return assembly_map
#mesh refinement and make h
def refine_mesh(assembly_map, n):
    cell_map = np.repeat(assembly_map, n, axis=0)
    cell_map = np.repeat(cell_map, n, axis = 1)
    h = ASSEMBLY_PITCH / n
    return cell_map, h
