import numpy as np
import math


def coes2rv(coes: tuple[float, float, float, float, float, float], body_mu: float, deg=False):
    # convert 6 classical orbital elements to position and velocity vectors in 3 d
    a, e, i, raan, aop, ta = coes
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
    # v_perifocal =


def eccentric_anomaly(anomaly, e, method='newton', tolerance=1e-8):
    if method == 'newton':
        M = anomaly
        pass
    elif method == 'tae':
        ta = anomaly
        pass

    return "no nothing implemented"
