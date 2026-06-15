#!/usr/bin/env python3

from fractions import Fraction

import h5py
import matplotlib.pyplot as plt
import numpy as np


# should be consistent with parameters.jl
# the projector workflow is canonical, so mu enters T but does not tune filling here
Lx         = 6
Ly         = 6
U          = 4.0
mu         = 0.0
Nbin       = 20
Nsweep     = 100
hs_channel = ":SU2"

# h5 data location
data = "analysis/data.h5"

try:
    f = h5py.File(data, "r")
except FileNotFoundError:
    raise RuntimeError(f"Missing analyzed data file: {data}")

# Load analyzed observables from the HDF5 output.
with f:
    sign = f["sign"][:, 0]
    ke = f["KE"][:, 0]
    pe = f["PE"][:, 0]
    szk = f["SzSzk"][:, :, 0, 0]

# The first axis stores [mean, error].
sign_mean, sign_err = np.real(sign[0]), np.real(sign[1])
ke_mean, ke_err = np.real(ke[0]), np.real(ke[1])
pe_mean, pe_err = np.real(pe[0]), np.real(pe[1])

# Reshape S(k) back to the 2D momentum grid.
szk_mean = np.real(szk[0, :]).reshape(Lx, Ly)
szk_err = np.real(szk[1, :]).reshape(Lx, Ly)

# The AFM ordering vector is (pi, pi) on the 6x6 grid.
afm_ix = Lx // 2
afm_iy = Ly // 2
peak_idx = np.unravel_index(np.argmax(szk_mean), szk_mean.shape)

# print basic informations
print(f"Read {data}")
print(f"L = {Lx}x{Ly}, U = {U}, mu = {mu}, Nbin = {Nbin}, Nsweep = {Nsweep}, hs_channel = {hs_channel}")
print(f"sign = {sign_mean:.6f} ± {sign_err:.6f}")
print(f"KE   = {ke_mean:.6f} ± {ke_err:.6f}")
print(f"PE   = {pe_mean:.6f} ± {pe_err:.6f}")
print(f"S(pi, pi) = {szk_mean[afm_ix, afm_iy]:.6f} ± {szk_err[afm_ix, afm_iy]:.6f}")
print(
    f"Peak S(k) at grid index ({peak_idx[0] + 1}, {peak_idx[1] + 1}) = "
    f"{szk_mean[peak_idx]:.6f} ± {szk_err[peak_idx]:.6f}"
)

# Plot the AFM structure factor heatmap in momentum space.
plt.figure(figsize=(6, 5))
im = plt.imshow(szk_mean.T, origin="lower", cmap="viridis", interpolation="nearest")
plt.colorbar(im, label=r"$S(\mathbf{k})$")
plt.scatter([afm_ix], [afm_iy], facecolors="none", edgecolors="white", s=120, linewidths=1.8)
tick_x = np.arange(Lx)
tick_y = np.arange(Ly)
tick_x_labels = [str(Fraction(2 * ix, Lx)) for ix in tick_x]
tick_y_labels = [str(Fraction(2 * iy, Ly)) for iy in tick_y]
plt.xticks(tick_x, tick_x_labels)
plt.yticks(tick_y, tick_y_labels)
plt.title(r"AFM structure factor in $k$ space")
plt.xlabel(r"$k_x / \pi$")
plt.ylabel(r"$k_y / \pi$")
plt.tight_layout()
plt.savefig("afm_structure_factor_heatmap.png", dpi=200)
plt.close()

print(f"Saved afm_structure_factor_heatmap.png")
