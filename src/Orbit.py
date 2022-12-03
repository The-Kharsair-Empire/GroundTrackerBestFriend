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

    E = eccentric_anomaly(ta, e, 'tae')

    r_norm = a * (1 - e ** 2) / (1 + e * math.cos(ta))

    r_perifocal = r_norm * np.array([math.cos(ta), math.sin(ta), 0])
    v_perifocal = math.sqrt(body_mu * a) / r_norm * np.array([-math.sin(E), math.cos(E) * math.sqrt(1 - e**2), 0])

    Cpe = perifocal2eci(raan, i, aop)

    r = np.dot(Cpe, r_perifocal)
    v = np.dot(Cpe, v_perifocal)

    return r, v


def eci2perifocal(raan, i, aop):
    s = math.sin
    c = math.cos
    return np.array([
        [-s(raan) * c(i) * s(aop) + c(raan) * c(aop), c(raan) * c(i) * s(aop) + s(raan) * c(aop), s(i) * s(aop)],
        [-s(raan) * c(i) * c(aop) - c(raan) * s(aop), c(raan) * c(i) * c(aop) - s(raan) * s(aop), s(i) * c(aop)],
        [s(raan) * s(i), -c(raan) * s(i), c(i)]
    ])


def perifocal2eci(raan, i, aop):
    return np.transpose(eci2perifocal(raan, i, aop))


def eccentric_anomaly(anomaly, e, method='newton', tolerance=1e-8):
    if method == 'newton':
        M = anomaly
        pass
    elif method == 'tae':
        ta = anomaly
        return 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(ta / 2.0))

    raise ArithmeticError("no eccentric anomaly implemented")
