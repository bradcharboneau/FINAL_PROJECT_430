from mesh import make_square_mesh, refine_mesh, ASSEMBLY_PITCH
from solver import build_matrices, power_iteration
from exact_answer import k_exact
from material_props import materials, FUEL
import numpy as np

n_assemblies = 8
refinement = [1,2,4,8,16,32,64]
L = n_assemblies * ASSEMBLY_PITCH
k_ex = k_exact(materials[FUEL], L)
h_vals = []
k_vals = []
err_vals = []
print(k_ex)
for n in refinement:
    cell_map, h = refine_mesh(make_square_mesh(n_assemblies), n)
    A1, A2, F1, F2, S12 = build_matrices(cell_map, h, materials)
    k, phi1, phi2, n_it = power_iteration(A1, A2, F1, F2, S12)
    err = abs(k - k_ex)
    h_vals.append(h)
    k_vals.append(k)
    err_vals.append(err)
    print(n, h, k, err, n_it)
