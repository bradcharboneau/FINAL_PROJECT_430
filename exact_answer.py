## exact answer
import numpy as np
from material_props import materials, FUEL
from mesh import ASSEMBLY_PITCH


def k_exact(xs, L):
    B2 = 2 * (np.pi / L)**2
    D1, D2 = xs["D"]
    sig_a1, sig_a2 = xs["sig_a"]
    nu_sig_f1, nu_sig_f2 = xs["nu_sig_f"]
    sig12 = xs["sig_s12"]
    R1 = sig_a1 + sig12 + D1 * B2
    R2 = sig_a2 + D2 * B2
    k = (nu_sig_f1 / R1) + (nu_sig_f2 * sig12) / (R1 * R2)
    return k
#test block
if __name__ == "__main__":
    xs = materials[FUEL]
    L = 3 * ASSEMBLY_PITCH
    print(k_exact(xs, L))
