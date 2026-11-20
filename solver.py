import numpy as np
import scipy.sparse as sp
#2d mapping start top left, cross row first
def cell_index(i, j, ny):
    return i* ny + j

#build out matrices for each cell

def build_matrices(cell_map, h, materials):
    n_rows, n_cols = cell_map.shape
    N = n_rows * n_cols
    A1 = sp.lil_matrix((N, N))
    A2 = sp.lil_matrix((N, N))
    F1 = np.zeros(N)
    F2 = np.zeros(N)
    S12 = np.zeros(N)
    for i in range(n_rows):
        for j in range(n_cols):
            p = cell_index(i, j, n_cols)
            mat_id = cell_map[i,j]
            xs = materials[mat_id]
            A1[p,p] += xs["sig_a"][0] + xs["sig_s12"]
            A2[p,p] += xs["sig_a"][1]
            F1[p] = xs["nu_sig_f"][0]
            F2[p] = xs["nu_sig_f"][1]
            S12[p] = xs["sig_s12"]
    return A1, A2, F1, F2, S12


##test block
if __name__ == "__main__":
    from mesh import make_square_mesh, refine_mesh
    from material_props import materials

    cell_map = make_square_mesh(2)
    A1, A2, F1, F2, S12 = build_matrices(cell_map, 1.0, materials)

    print(A1.shape, F1.shape)
    print(A1.toarray())
    print(F2)
