import numpy as np
import math


# TODO: try other different lambert solver algorithms: Gauss's method
def lambert_solver_pykey_project(r0: np.ndarray, r1: np.ndarray, transfer_time, mu, prograde=True):
    # TODO: I mean this is totally wrong, no idea
    r0_norm = np.linalg.norm(r0)
    r1_norm = np.linalg.norm(r1)
    c = np.linalg.norm((r0 - r1))
    s = (r0_norm + r1_norm + c) / 2
    lambda_ = math.sqrt(1 - c / s)
    t = transfer_time * math.sqrt(2 * mu / s ** 3)

    ir0 = r0 / r0_norm
    ir1 = r1 / r1_norm

    ih = np.cross(ir0, ir1)
    ih = ih / np.linalg.norm(ih)

    it0 = np.cross(ih, ir0)
    it0 = it0 / np.linalg.norm(it0)

    it1 = np.cross(ih, ir1)
    it1 = it1 / np.linalg.norm(it1)

    if ih[1] < 0 != prograde:
        it0 = -it0
        it1 = -it1
        lambda_ = -lambda_

    x = iterative_root_finder(lambda_, t)

    y = math.sqrt(1 - lambda_ ** 2 * (1 - x ** 2))
    z = lambda_ * y
    rho = (r0_norm - r1_norm) / c
    gamma = math.sqrt(mu * s / 2)

    vr0 = (z - x) - rho * (x + z)
    vr1 = (x - z) - rho * (x + z)
    vt = math.sqrt(1 - rho ** 2) * (y + lambda_ * x)

    v0 = (gamma / r0_norm) * (vr0 * ir0 + vt * it0)
    v1 = (gamma / r1_norm) * (vr1 * ir1 + vt * it1)

    return v0, v1


def iterative_root_finder(lambda_, t):
    x = initial_guess(lambda_, t)
    delta = 1
    iterations = 0
    while abs(delta) >= 0.00001 or iterations == 15:
        delta = householders_method(lambda_, t, x)
        x -= delta
        iterations += 1

    return x


def initial_guess(lambda_, t):

    d2r = np.pi / 180.0
    t0 = abs(d2r * math.acos(lambda_) + lambda_ * math.sqrt(1 - lambda_ ** 2))
    # adding a abs() is the only way to fix it, but it seems that the direction is still reverse
    t1 = (2 / 3) * (1 - lambda_ ** 3)
    # print(t, t0, t1)

    if t >= t0:
        return (t0 / t) ** (2 / 3) - 1     # warning: risk of complex value
    elif t <= t1:
        # print(2)
        # print((5 * t1 * (t1 - t)) / (2 * t * (1 - lambda_ ** 5)) + 1)
        return (5 * t1 * (t1 - t)) / (2 * t * (1 - lambda_ ** 5)) + 1
    else:
        # print(3)
        # print((t0 / t) ** (math.log(t1 / t0) / math.log(2)) - 1)
        return (t0 / t) ** (math.log(t1 / t0) / math.log(2)) - 1


def householders_method(lambda_, t, x):
    a = 1 - x ** 2

    y = math.sqrt(1 - lambda_ ** 2 * a)
    tau = time_of_flight(lambda_, a, x, y)
    delta = tau - t

    dt = (3 * tau * x - 2 + 2 * (lambda_ ** 3) * x / y) / a
    ddt = (3 * tau + 5 * x * dt + 2 * (1 - lambda_ ** 2) * (lambda_ ** 3) / (y ** 3)) / a
    dddt = (7 * x * ddt + 8 * dt - 6 * (1 - lambda_ ** 2) * (lambda_ ** 5) * x / (y ** 5)) / a

    return delta * (dt ** 2 - delta * ddt / 2) / (dt * (dt ** 2 - delta * ddt) + (dddt * delta ** 2) / 6)


def time_of_flight(lambda_, a, x, y):
    d2r = np.pi / 180.0
    b = math.sqrt(abs(a))
    f = b * (y - lambda_ * x)
    g = lambda_ * a + x * y
    psi = d2r * math.acos(g) if a > 0 else math.log(max(1e-300, f + g), math.e)
    return (psi / b - x + lambda_ * y) / a


def lambert_universal_variable_algorithm(r0: np.ndarray, r1: np.ndarray,
                                         transfer_time, mu, prograde=True, tol=1e-6,
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
        if abs(transfer_time - dt_tilda) < tol:
            solved = True
            break

        if dt_tilda <= transfer_time:
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
