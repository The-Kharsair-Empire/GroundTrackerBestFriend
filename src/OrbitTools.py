import numpy as np
import math


def coes2rv(a, e, i, raan, aop, ta, body_mu: float, deg=False):
    # convert 6 classical orbital elements to position and velocity vectors in 3 d
    # semi-major Axis,
    # Eccentricity,
    # Inclination,
    # Right Ascension of the Ascending Node
    # Argument Of Periapsis,
    # True Anomaly
    if deg:
        d2r = np.pi / 180.0
        i *= d2r
        raan *= d2r
        aop *= d2r
        ta *= d2r

    E = ta2E(ta, e)

    r_norm = a * (1 - e ** 2) / (1 + e * math.cos(ta))

    r_perifocal = r_norm * np.array([math.cos(ta), math.sin(ta), 0])
    v_perifocal = math.sqrt(body_mu * a) / r_norm * np.array([-math.sin(E), math.cos(E) * math.sqrt(1 - e ** 2), 0])

    Cpe = perifocal2eci(raan, i, aop)

    r = np.dot(Cpe, r_perifocal)
    v = np.dot(Cpe, v_perifocal)

    return r, v


def rv2coes(r, v, body_mu, ta_in_time=False, t=None, deg=False):

    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)

    # a = -body_mu / 2 * (v_norm ** 2 / 2 - body_mu / r_norm)

    e_vector = np.cross(v, h) / body_mu - r / r_norm
    e_vector_2 = ((np.linalg.norm(v) ** 2 - body_mu / r_norm) * r - np.dot(r, v) * v) / body_mu
    e = np.linalg.norm(e_vector)
    e_2 = np.linalg.norm(e_vector_2)

    i = np.arccos(h[2] / h_norm)

    # n vector is pointing from center of body to equatorial ascending node
    n_vector = np.cross([0, 0, 1], h)
    n_norm = np.linalg.norm(n_vector)

    if np.isclose(n_norm, 0):
        print("something wrong with n_norm, too close to 0")
        print("because h is too close to z [0, 0, 1], the angle between them almost 0, orbital plane too equatorial")
        print(f"angle between them: {np.arccos(np.dot([0, 0, 1], h) / (np.linalg.norm([0, 0, 1]) * h_norm))}")
        print(f"checking n_vector: {n_vector}")
        print(f"checking h: {h}")
        print(f"checking r: {r}")
        print(f"checking v: {v}")
        n_norm = 0.01

    raan = np.arccos(n_vector[0] / n_norm)  # n_vector[0] = np.dot([1, 0, 0], n_vector)
    if n_vector[1] < 0:
        raan = 2 * np.pi - raan

    aop = np.arccos(np.dot(n_vector, e_vector) / (n_norm * e))
    if e_vector[2] < 0:
        aop = 2 * np.pi - aop

    ta = np.arccos(np.dot(e_vector, r) / (e * r_norm))
    if np.dot(r, v) < 0:
        ta = 2 * np.pi - ta

    a = r_norm * (1 + e * np.cos(ta)) / (1 - e ** 2)
    # a_2 = r_norm * (1 + e * np.cos(ta)) / (1 - e ** 2)

    # if not np.isclose(a, a_2):
    #     print(f"two algorithm for semi-major axis doesn't match: {a} {a_2}")

    if not np.isclose(e, e_2):
        print(f"two algorithm for eccentricity doesn't match: {e} {e_2}")

    tp = None
    if ta_in_time:
        if t is None:
            raise ArithmeticError("Should provide time if you want the result with tp")
        E = ta2E(ta, e)
        tp = t - (E - e * np.sin(E)) / np.sqrt(body_mu / a ** 3)

    if deg:
        r2d = 180.0 / np.pi
        i *= r2d
        raan *= r2d
        aop *= r2d
        ta *= r2d

    if ta_in_time:
        return a, e, i, raan, aop, tp
    else:
        return a, e, i, raan, aop, ta


def eci2perifocal(raan, i, aop):
    s = np.sin
    c = np.cos
    return np.array([
        [-s(raan) * c(i) * s(aop) + c(raan) * c(aop), c(raan) * c(i) * s(aop) + s(raan) * c(aop), s(i) * s(aop)],
        [-s(raan) * c(i) * c(aop) - c(raan) * s(aop), c(raan) * c(i) * c(aop) - s(raan) * s(aop), s(i) * c(aop)],
        [s(raan) * s(i), -c(raan) * s(i), c(i)]
    ])


def perifocal2eci(raan, i, aop):
    return np.transpose(eci2perifocal(raan, i, aop))


def ta2E(ta, e):
    # getting eccentric anomaly from true anomaly
    return 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(ta / 2.0))


def M2E(M, e, tolerance=1e-8, max_step=200):
    # newton's method for solving eccentric anomaly using mean anomaly
    def iterative_step(E):
        return (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

    if M < np.pi / 2.0:
        E0 = M + e / 2.0
    else:
        E0 = M - e

    E1 = E0 - iterative_step(E0)
    for n in range(200):
        if abs(E1 - E0) < tolerance:
            return E1
        else:
            E0 = E1
            E1 -= iterative_step(E1)


def E2ta(E, e):
    # getting true anomaly from eccentric anomaly
    return 2 * np.arctan(math.sqrt((1 + e) / (1 - e)) * np.tan(E / 2.0))


def get_r_mag(a, e, E):
    # position from semi-major axis, eccentricity and Eccentric anomaly
    return a * (1 - e * np.cos(E))
