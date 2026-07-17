#!/usr/bin/env python3
"""3-D wealth-age density surfaces: status quo vs the largest-grant steady state.

Produces two figures (linear pdf and log10 pdf), each with two smooth 3-D
surfaces of the joint density f(a, h) = sum_y m / (w_a * dh), under step2
(Destatis) mortality:
    results/step2/density3d_pdf.png
    results/step2/density3d_logpdf.png

Input preference: results/step2/maxgrant_step2.mat (true largest balanceable
grant, from run_max_grant on the remote). Falls back to the G=200k steady
state (results/step2/flatgrant_200k_step2.mat, ~94% of the frontier) with a
clear label if the max-grant run is not available yet.

Run:  python3 postprocess/make_density3d.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


def grids_from_params(p):
    Na, Ny = int(p.Na), int(p.Ny)
    a = float(p.a_max) * np.linspace(0, 1, Na) ** float(p.a_curv)
    h = np.arange(int(p.h0), int(p.h_max) + 1)
    Nh = h.size
    da = np.diff(a)
    wa = np.empty(Na)
    wa[0], wa[-1] = da[0] / 2, da[-1] / 2
    wa[1:-1] = (da[:-1] + da[1:]) / 2
    return a, h, wa, Na, Ny, Nh


def density_ah(m, p):
    """Joint density f(a, h): masses summed over y, divided by cell measure."""
    a, h, wa, Na, Ny, Nh = grids_from_params(p)
    M = m.reshape((Na, Ny, Nh), order="F")     # column-major (MATLAB) layout
    f = M.sum(axis=1) / wa[:, None]            # dh = 1
    return a, h, f


def load_solutions():
    sq_file = os.path.join(BASE, "step2", "flatgrant_20k_step2.mat")
    S = loadmat(sq_file, squeeze_me=True, struct_as_record=False)
    p_sq, m_sq = S["out"].sol_sq.params, S["out"].sol_sq.m

    mg_file = os.path.join(BASE, "step2", "maxgrant_step2.mat")
    if os.path.exists(mg_file):
        Smg = loadmat(mg_file, squeeze_me=True, struct_as_record=False)
        p_g, m_g = Smg["params"], Smg["m_grant"]
        label = f"max grant ({Smg['out'].G_max_eur/1e3:.0f}k EUR @ surtax cap)"
    else:
        Sg = loadmat(os.path.join(BASE, "step2", "flatgrant_200k_step2.mat"),
                     squeeze_me=True, struct_as_record=False)
        p_g, m_g = Sg["out"].sol.params, Sg["out"].sol.m
        label = "Grunderbe 200k EUR (near-frontier stand-in)"
    return (p_sq, m_sq), (p_g, m_g), label


LOG_FLOOR = 1e-10   # density floor for the log plot


def surface(ax, a, h, f, eur, a_cap_keur, log=False, title="", zlab="", zclip=None):
    a_keur = a * eur / 1e3
    keep = a_keur <= a_cap_keur
    A, H = np.meshgrid(a_keur[keep], h, indexing="ij")
    Z = f[keep, :].copy()
    if log:
        Z = np.log10(np.maximum(Z, LOG_FLOOR))
    elif zclip is not None:
        Z = np.minimum(Z, zclip)
    ax.plot_surface(H, A, Z, cmap="viridis", rstride=1, cstride=1,
                    linewidth=0, antialiased=True)
    ax.set_xlabel("age")
    ax.set_ylabel("wealth (1000 EUR)")
    ax.set_zlabel(zlab, labelpad=8)
    ax.set_title(title, pad=10)
    ax.view_init(elev=28, azim=-135)


def main():
    (p_sq, m_sq), (p_g, m_g), glabel = load_solutions()
    eur = float(p_sq.eur_per_unit)
    a, h, f_sq = density_ah(m_sq, p_sq)
    _, _, f_g = density_ah(m_g, p_g)

    # clip the linear pdf at the largest density found ABOVE 100k EUR wealth,
    # so the old-age near-zero-wealth spike does not flatten all structure
    hi_region = a * eur / 1e3 >= 100
    zclip = 1.05 * max(f_sq[hi_region, :].max(), f_g[hi_region, :].max())

    for log, cap_keur, fname in [(False, 1500, "density3d_pdf.png"),
                                 (True, 5000, "density3d_logpdf.png")]:
        fig = plt.figure(figsize=(14, 6))
        zl = "log10 density" if log else "density f(a,h)"
        ax1 = fig.add_subplot(1, 2, 1, projection="3d")
        surface(ax1, a, h, f_sq, eur, cap_keur, log, "status quo", zl, zclip)
        ax2 = fig.add_subplot(1, 2, 2, projection="3d")
        surface(ax2, a, h, f_g, eur, cap_keur, log, glabel, zl, zclip)
        # common z-limits for comparability
        lo = min(ax1.get_zlim()[0], ax2.get_zlim()[0])
        hi = max(ax1.get_zlim()[1], ax2.get_zlim()[1])
        ax1.set_zlim(lo, hi); ax2.set_zlim(lo, hi)
        note = " — log scale" if log else f" — z clipped at {zclip:.2g} (old-age a≈0 spike)"
        fig.suptitle("Wealth-age density, step 2 (Destatis) mortality" + note)
        out = os.path.join(BASE, "step2", fname)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("saved", out)


if __name__ == "__main__":
    main()
