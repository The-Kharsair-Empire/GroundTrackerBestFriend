from dataclasses import dataclass
from Utility import NoneRefersDefault, DefaultVal
import numpy as np


@dataclass
class CelestialBody(NoneRefersDefault):
    name: str
    mass: float  # kg
    mu: float  # km^3 / s^2
    radius: float  # km
    j2: float = DefaultVal(0.0)  # Legendre polynomial -C(2, 2)
    angular_velocity: np.ndarray = DefaultVal(np.zeros(3))  # rad / s
    rhos: np.ndarray = DefaultVal(np.zeros(3))
    zs: np.ndarray = DefaultVal(np.zeros(3))


sun = CelestialBody('Sun', 1.989e30, 1.32712e11, 695700.0)

earth = CelestialBody('Earth', 5.972e24, 3.986e5, 6378.0, 0.001082635854
                      , np.array([0.0, 0.0, 72.9211e-6])
                      , np.array([2.059e-4, 5.909e-11, 3.561e-15]) * 10 ** 8
                      , np.array([63.096, 251.189, 1000.0]))
