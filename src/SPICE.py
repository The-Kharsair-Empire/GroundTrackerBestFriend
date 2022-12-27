import spiceypy as spice
import numpy as np
import os

loaded_spice_kernel = []
default_lsk_kernel = "latest_leapseconds.tls.pc"


def spice_load_kernels(kernel_files: list[str]):
    assert kernel_files is not None, "Kernel file is None"

    if isinstance(kernel_files, str):
        kernel_files = [kernel_files]

    assert len(kernel_files) > 0 and isinstance(kernel_files[0], str), "Kernel files invalid"

    kernel_files.append(default_lsk_kernel)

    for kernel in kernel_files:
        if kernel not in loaded_spice_kernel:
            spice.furnsh(os.path.join(os.path.dirname(__file__), '..', 'file', 'spice_data', kernel))
            loaded_spice_kernel.append(kernel)
        else:
            print(kernel, "is loaded, skipping")


def spice_get_spk_objects(filename, verbose=False):
    filename = os.path.join(os.path.dirname(__file__), '..', 'file', 'spice_data', filename)
    objects = spice.spkobj(filename)
    ids, names, tcs_sec, tcs_cal = [], [], [], []

    n = 0

    if verbose:
        print(f'\nObjects in {filename}')

    for o in objects:
        ids.append(o)

        tc_sec = spice.wnfetd(spice.spkcov(filename, ids[n]), n)
        tc_format = 'DD/MON/YYYY HR:MN:SC.### (TDB) ::TDB'
        tc_cal = (spice.timout(tc_sec[0], tc_format), spice.timout(tc_sec[1], tc_format))
        # tc_cal = [spice.timout(f, 'DD/MON/YYYY HR:MN:SC.### (TDB) ::TDB') for f in tc_sec]

        tcs_sec.append(tc_sec)
        tcs_cal.append(tc_cal)

        try:
            names.append(spice.bodc2n(o))  # return name of body given spice id
        except:
            names.append('Unknown Object')  # some random body might not have name in spice database

        if verbose:
            print(f'id: {ids[-1]}\t\tname: {names[-1]}\t\ttc: {tc_cal[0]} --> {tc_cal[1]}')

    return ids, names, tcs_sec, tcs_cal


def spice_get_ephemeris_data(target, times, frame, observer, abcorr='NONE'):
    return np.array(spice.spkezr(target, times, frame, abcorr, observer)[0])


def tc2array(tcs, steps):
    arr = np.zeros((steps, 1))
    arr[:, 0] = np.linspace(tcs[0], tcs[1], steps)
    return arr


def rv2coes(state, et, mu, deg=False):
    rp, e, i, raan, aop, ma, t0, mu, ta, a, T = spice.oscltx(state, et, mu)

    if deg:
        r2d = 180.0 / np.pi
        i *= r2d
        ta *= r2d
        aop *= r2d
        raan *= r2d

    return a, e, i, raan, aop, ta

