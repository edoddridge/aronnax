"""Regenerate the figures in docs/figs used by verification.rst.

Run from the repository root:  python docs/make_verification_figures.py
"""
import os.path as p
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = p.dirname(p.abspath(__file__))
sys.path.insert(0, p.join(here, "..", "test"))
import physics_test as pt

figs = p.join(here, "figs")
DXS = [20e3, 10e3, 5e3, 2.5e3]
F = 1e-4


def run(function, *args, **kwargs):
    with tempfile.TemporaryDirectory() as d:
        return function(d, *args, **kwargs)


def error_by_resolution():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
    for ax, f in zip(axes, [0., F]):
        for layers in [1, 2, 3, 4]:
            errors = [pt.relative_error(*run(pt.run_wind_forcing, layers, dx, f)[1:])
                      for dx in DXS]
            ax.loglog(np.array(DXS) / 1e3, errors, "o-",
                      label="n-layer, {} layer{}".format(
                          layers, "s" if layers > 1 else ""))
        errors = [pt.relative_error(
            *run(pt.run_wind_forcing, 1, dx, f, red_grav=True)[1:])
            for dx in DXS]
        ax.loglog(np.array(DXS) / 1e3, errors, "k--", label="reduced gravity")
        if f:
            ax.set_yscale("linear")
            ax.set_ylim(2e-5, 3e-5)
        ax.set_xlabel("Grid spacing (km)")
        ax.set_title("f = {:g} 1/s".format(f))
    axes[0].set_ylabel("Maximum relative momentum error")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(p.join(figs, "verification_wind_error.png"), dpi=150)
    plt.close(fig)


def momentum_evolution():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (layers, red_grav) in zip(axes, [(4, False), (1, True)]):
        t, sim, exact = run(pt.run_wind_forcing, layers, 10e3, F,
                            red_grav=red_grav, n_dumps=50)
        t = t / 86400.
        ax.plot(t, exact.real, "-", label="Expected, zonal")
        ax.plot(t, exact.imag, "-", label="Expected, meridional")
        ax.plot(t, sim.real, "o", ms=4, mfc="none", label="Simulated, zonal")
        ax.plot(t, sim.imag, "s", ms=4, mfc="none", label="Simulated, meridional")
        ax.set_xlabel("Time (days)")
        ax.set_title("reduced gravity" if red_grav
                     else "n-layer, {} layers".format(layers))
    axes[0].set_ylabel("Momentum")
    axes[0].set_ylim(-4.2e7, 2.8e7)
    axes[0].legend(loc="lower left", ncol=2, fontsize=7)
    fig.tight_layout()
    fig.savefig(p.join(figs, "verification_wind_momentum.png"), dpi=150)
    plt.close(fig)


def conservation():
    fig, ax = plt.subplots(figsize=(6, 4))
    for layers, red_grav in [(1, False), (2, False), (3, False), (4, False),
                             (4, True)]:
        its, mom = run(pt.run_unforced, layers, red_grav)
        change = np.abs(mom - mom[0]) / np.abs(mom[0])
        ax.semilogy(its[1:], change[1:], "o-",
                    label="{}, {} layer{}".format(
                        "reduced gravity" if red_grav else "n-layer", layers,
                        "s" if layers > 1 else ""))
    ax.set_xlabel("Time step")
    ax.set_ylabel("Relative change in momentum")
    ax.legend()
    fig.tight_layout()
    fig.savefig(p.join(figs, "verification_conservation.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    error_by_resolution()
    momentum_evolution()
    conservation()
