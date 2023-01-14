import numpy as np
import math

# TODO: try other different lambert solver algorithms: Gauss's method, PyKep project


def lambert_universal_variable_algorithm(r0: np.ndarray, r1: np.ndarray,
                                         dt, mu, prograde=True, tol=1e-6,
                                         max_steps=200,
                                         psi_0=0, psi_u=4 * np.pi ** 2, psi_l=-4 * np.pi):
    if prograde:
        tm = 1
    else:
        tm = -1
    sqrt_mu = np.sqrt(mu)
    r0_norm = np.linalg.norm(r0)
    r1_norm = np.linalg.norm(r1)
    gamma = np.dot(r0, r1) / (r0_norm * r1_norm)

    beta = tm * np.sqrt(1 - gamma ** 2)
    A = tm * np.sqrt(r0_norm * r1_norm * (1 + gamma))
    if A == 0:
        print("cannot calculate solution")
        return np.array([0, 0, 0]), np.array([0, 0, 0])

    step = 0
    B = 0
    solved = False

    psi = psi_0

    c2 = 0.5
    c3 = 1 / 6.0

    for n in range(max_steps):

        B = r0_norm + r1_norm + A * (psi * c3 - 1) / np.sqrt(c2)

        if A > 0.0 > B:
            psi_l += np.pi
            B *= -1

        chi3 = np.sqrt(B / c2) ** 3
        dt_tilda = (chi3 * c3 + A * np.sqrt(B)) / sqrt_mu
        if abs(dt - dt_tilda) < tol:
            solved = True
            break

        if dt_tilda <= dt:
            psi_l = psi
        else:
            psi_u = psi

        psi = (psi_u + psi_l) / 2.0
        c2 = get_c2(psi)
        c3 = get_c3(psi)

    if not solved:
        print('algorithm did not converge')
        return np.array([0, 0, 0]), np.array([0, 0, 0])

    F = 1 - B / r0_norm
    G = A * np.sqrt(B / mu)
    Gdot = 1 - B / r1_norm

    v0 = (r1 - r0 * F) / G
    v1 = (Gdot * r1 - r0) / G

    return v0, v1


def get_c2(psi):
    return (1 - math.cos(math.sqrt(psi))) / psi


def get_c3(psi):
    return (math.sqrt(psi) - math.sin(math.sqrt(psi))) / (psi * math.sqrt(psi))


def get_c2psi_c3psi(psi, psi_m):
    K = 0
    while not abs(psi) <= psi_m:
        K += 1
        psi /= 4

    c2, c3 = algorithm_1()
    while K != 0:
        K -= 1
        c2, c3 = algorithm_other()
        psi *= 4

    return c2, c3


def algorithm_1():
    return 0, 0


def algorithm_other():
    return 0, 0
