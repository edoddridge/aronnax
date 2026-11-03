"""Momentum verification for the n-layer and reduced gravity physics.

A localised wind stress is applied to the top layer of a doubly periodic
f-plane with all explicit viscosity and drag switched off. Pressure forces
and the exchange of momentum between layers integrate to zero over a
doubly periodic domain, so the domain-integrated, signed momentum summed
over all layers obeys

    dM/dt = F + f k x M,

where F is the total wind force. With Z = Mx + i My this is
dZ/dt = F (1 + i) - i f Z, which has a closed-form solution.

The signed sum is essential. The old verification script summed
absolute values, which also counts the momentum of radiated gravity waves
and exaggerates the n-layer result (issue #100).
"""

import glob
import os.path as p

import numpy as np
import pytest

import aronnax as aro
import aronnax.driver as drv

RHO0 = 1035.
WIND = 1e-5
NX = 20
DT = 100.
N_TIME_STEPS = 1000

CONFIG_TEMPLATE = """
[numerics]
au = 0.
kh = 0.
ar = 0.
bot_drag = 0.
dt = {dt}
slip = 0.
n_time_steps = {n_time_steps}
dump_freq = {dump_freq}
av_freq = 0.
hmin = 100
maxits = 1000
eps = 1e-9
freesurf_fac = 1.
thickness_error = 1e-2

[model]
hmean = {hmean}
h0 = {h0}
red_grav = {red_grav}

[pressure_solver]
n_proc_x = 1
n_proc_y = 1

[physics]
g_vec = {g_vec}
rho0 = 1035.

[grid]
nx = {nx}
ny = {nx}
layers = {layers}
dx = {dx}
dy = {dx}
wet_mask_file = :rectangular_pool:

[external_forcing]
dump_wind = no
relative_wind = no
"""


def write_config(work_dir, layers, dx, red_grav):
    thickness = 400.
    g_vec = ", ".join(["9.8"] + ["0.01"] * (layers - 1))
    if red_grav:
        g_vec = ", ".join(["0.01"] * layers)
    with open(p.join(str(work_dir), "aronnax.conf"), "w") as conf:
        conf.write(CONFIG_TEMPLATE.format(
            dt=DT, n_time_steps=N_TIME_STEPS, dump_freq=DT * N_TIME_STEPS / 5,
            hmean=thickness * layers, h0=thickness * layers,
            red_grav="yes" if red_grav else "no", g_vec=g_vec, nx=NX,
            layers=layers, dx=dx))
    return thickness


def wet_mask(X, Y):
    return np.ones(X.shape)


def read_momentum(work_dir, hfile, layers, dx):
    """Signed momentum (Mx, My) summed over all layers."""
    h = aro.interpret_raw_file(hfile, NX, NX, layers)
    u = aro.interpret_raw_file(hfile.replace("snap.h", "snap.u"),
                               NX, NX, layers)
    v = aro.interpret_raw_file(hfile.replace("snap.h", "snap.v"),
                               NX, NX, layers)
    mx = dx * dx * RHO0 * np.sum(h * (u[:, :, 1:] + u[:, :, :-1]) / 2.)
    my = dx * dx * RHO0 * np.sum(h * (v[:, 1:, :] + v[:, :-1, :]) / 2.)
    return mx + 1j * my


def run_wind_forcing(work_dir, layers, dx, f, red_grav=False):
    """Run the forced simulation and return times, simulated and exact Z."""
    nx = NX
    thickness = write_config(work_dir, layers, dx, red_grav)

    def wind(X, Y):
        w = np.zeros(Y.shape)
        w[nx // 2, nx // 2] = WIND
        return w

    drv.simulate(work_dir=str(work_dir), exe="aronnax_test",
                 init_h_file=[thickness] * layers,
                 zonal_wind_file=[wind], meridional_wind_file=[wind],
                 wet_mask_file=[wet_mask],
                 f_u_file=[f], f_v_file=[f],
                 nx=nx, ny=nx, dx=dx, dy=dx)

    force = WIND * dx * dx
    times, simulated, exact = [], [], []
    for hfile in sorted(glob.glob(p.join(str(work_dir), "output/snap.h.*"))):
        iteration = int(hfile[-10:])
        simulated_z = read_momentum(work_dir, hfile, layers, dx)
        # The first Adams-Bashforth steps are forward Euler, which shifts
        # the effective time by two steps.
        t = (iteration + 2) * DT
        if f == 0.:
            z = force * (1 + 1j) * t
        else:
            z = force * (1 + 1j) * (1 - np.exp(-1j * f * t)) / (1j * f)
        times.append(t)
        simulated.append(simulated_z)
        exact.append(z)
    return np.array(times), np.array(simulated), np.array(exact)


def relative_error(simulated, exact):
    return np.max(np.abs(simulated - exact)) / np.max(np.abs(exact))


@pytest.mark.parametrize("layers", [1, 2, 3, 4])
@pytest.mark.parametrize("dx", [20e3, 10e3, 5e3])
def test_wind_momentum_without_rotation(tmp_path, layers, dx):
    _, simulated, exact = run_wind_forcing(tmp_path, layers, dx, f=0.)
    assert len(simulated) > 3
    assert relative_error(simulated, exact) < 1e-7


@pytest.mark.parametrize("layers", [1, 2, 3, 4])
@pytest.mark.parametrize("dx", [20e3, 10e3, 5e3])
def test_wind_momentum_on_f_plane(tmp_path, layers, dx):
    _, simulated, exact = run_wind_forcing(tmp_path, layers, dx, f=1e-4)
    assert relative_error(simulated, exact) < 1e-3


def test_wind_momentum_error_independent_of_resolution(tmp_path):
    errors = []
    for k, dx in enumerate([20e3, 10e3, 5e3]):
        work_dir = tmp_path / str(k)
        work_dir.mkdir()
        _, simulated, exact = run_wind_forcing(work_dir, 2, dx, f=1e-4)
        errors.append(relative_error(simulated, exact))
    assert max(errors) < 1e-3
    assert max(errors) - min(errors) < 1e-4


def test_n_layer_matches_reduced_gravity_momentum(tmp_path):
    n_dir = tmp_path / "n"
    r_dir = tmp_path / "r"
    n_dir.mkdir()
    r_dir.mkdir()
    _, n_layer, _ = run_wind_forcing(n_dir, 1, 10e3, f=1e-4)
    _, red_grav, _ = run_wind_forcing(r_dir, 1, 10e3, f=1e-4, red_grav=True)
    assert np.max(np.abs(n_layer - red_grav)) / np.max(np.abs(red_grav)) < 1e-3


@pytest.mark.parametrize("layers", [1, 2, 3, 4])
@pytest.mark.parametrize("red_grav", [False, True])
def test_unforced_momentum_is_conserved(tmp_path, layers, red_grav):
    if red_grav and layers == 1:
        pytest.skip("covered by the n-layer case")
    dx = 10e3
    thickness = write_config(tmp_path, layers, dx, red_grav)

    def init_u(X, Y):
        return 3e-5 * np.exp(-((X - X.mean()) ** 2 + (Y - Y.mean()) ** 2)
                             / (2 * (3 * dx) ** 2))

    drv.simulate(work_dir=str(tmp_path), exe="aronnax_test",
                 init_h_file=[thickness] * layers,
                 init_u_file=[init_u] + [0.] * (layers - 1),
                 wet_mask_file=[wet_mask],
                 f_u_file=[0.], f_v_file=[0.],
                 nx=NX, ny=NX, dx=dx, dy=dx)
    hfiles = sorted(glob.glob(p.join(str(tmp_path), "output/snap.h.*")))
    momentum = np.array([read_momentum(tmp_path, hf, layers, dx)
                         for hf in hfiles])
    assert len(momentum) > 3
    assert np.max(np.abs(momentum - momentum[0])) / np.abs(momentum[0]) < 1e-10
