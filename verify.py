from mesh import make_square_mesh, refine_mesh, ASSEMBLY_PITCH
from solver import build_matrices, power_iteration
from exact_answer import k_exact
from material_props import materials, FUEL
import numpy as np
import matplotlib.pyplot as plt

n_assemblies = 8
refinement = [1,2,4,8,16,32]
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
orders = []
for i in range(len(err_vals) - 1):
    order = np.log(err_vals[i] / err_vals[i+1]) / np.log(h_vals[i] / h_vals[i+1])
    orders.append(order)
    print(f"h = {h_vals[i]:.5f} to {h_vals[i+1]:.5f}:  order = {order:.5f}")

fig, ax = plt.subplots()
h_arr = np.array(h_vals)
ref = err_vals[-1] * (h_arr / h_arr[-1])**2
#fitting line for reference
ax.loglog(h_arr, err_vals, "o-", label="FV solver")
ax.loglog(h_arr, ref, "--", color="red", label="$h^2$ reference")

ax.set_xlabel("cell size h (cm)")
ax.set_ylabel(r"$abs(k-k_{exact})$")
ax.set_title(f"Bare square, L = {L:.0f} cm")
ax.grid(True, which="both", alpha = 0.3)
ax.legend(loc="upper left")
ax.text(0.95, 0.05, f"observed order p = {orders[-1]:.5f}",
        transform=ax.transAxes, ha="right")

plt.show()
fig.tight_layout()
fig.savefig("verification_order.png", dpi=200)
