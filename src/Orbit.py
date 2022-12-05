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
    v_perifocal = math.sqrt(body_mu * a) / r_norm * np.array([-math.sin(E), math.cos(E) * math.sqrt(1 - e**2), 0])

    Cpe = perifocal2eci(raan, i, aop)

    r = np.dot(Cpe, r_perifocal)
    v = np.dot(Cpe, v_perifocal)

    return r, v


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
    return a * (1 - e * np.cos(E))
