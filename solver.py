import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

#2d mapping start top left, cross row first
def cell_index(i, j, ny):
    return i* ny + j
#harmonic average solver for face
def harmonic_mean(D_P, D_Q):
    return 2 * D_P * D_Q / (D_P + D_Q)

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
            for di, dj in [(-1,0), (1,0), (0,-1), (0,1)]:
                ni, nj = i + di, j + dj
                if 0 <= ni < n_rows and 0 <= nj < n_cols:
                    q = cell_index(ni, nj, n_cols)
                    xs_q = materials[cell_map[ni, nj]]
                    D_face = harmonic_mean(xs["D"][0], xs_q["D"][0])
                    A1[p,p] += D_face / h**2
                    A1[p,q] -= D_face / h**2
                    D_face = harmonic_mean(xs["D"][1], xs_q["D"][1])
                    A2[p,p] += D_face / h**2
                    A2[p,q] -= D_face / h**2
                #dirichlet boundary here for simplicity
                else:
                    A1[p,p] += 2 * xs["D"][0] / h**2
                    A2[p,p] += 2 * xs["D"][1] / h**2

    return A1.tocsc(), A2.tocsc(), F1, F2, S12
#power iteration time, arbitrary stopping criteria
def power_iteration(A1, A2, F1, F2, S12, tol_k=1e-8, tol_Q=1e-6, max_iter=1000):
    lu1 = spla.splu(A1)
    lu2 = spla.splu(A2)



##test block
if __name__ == "__main__":
    from mesh import make_square_mesh, refine_mesh
    from material_props import materials

    cell_map = make_square_mesh(2)
    A1, A2, F1, F2, S12 = build_matrices(cell_map, 1.0, materials)

    print(A1.shape, F1.shape)
    print(A1.toarray())
    print(F2)
