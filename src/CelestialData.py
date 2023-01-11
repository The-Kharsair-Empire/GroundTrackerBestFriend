from dataclasses import dataclass
import numpy as np
from dataclasses import dataclass, fields, field


d2r = np.pi / 180.0
r2d = 180.0 / np.pi
km2AU = 1.4959787e8
day2sec = 86400


@dataclass
class CelestialBody:
    name: str
    mass: float  # kg
    mu: float  # km^3 / s^2
    radius: float  # km
    j2: float = 0.0  # Legendre polynomial -C(2, 2)
    angular_velocity: np.ndarray = field(default_factory=lambda:np.zeros(3))  # rad / s
    rhos: np.ndarray = field(default_factory=lambda:np.zeros(3))
    zs: np.ndarray = field(default_factory=lambda:np.zeros(3))
    dist2parent: float = 0.0
    orbital_period: float = 0.0
    G1: float = 0.0  # kg-km^3

    def __post_init__(self):
        for field in fields(self):
            if getattr(self, field.name) is None:
                setattr(self, field.name, field.default)





sun = CelestialBody('Sun', 1.989e30, 1.32712e11, 695700.0, G1=1e8)

earth = CelestialBody('Earth', 5.972e24, 3.986e5, 6378.0, 0.001082635854
                      , np.array([0.0, 0.0, 72.9211e-6])
                      , np.array([2.059e-4, 5.909e-11, 3.561e-15]) * 10 ** 8
                      , np.array([63.096, 251.189, 1000.0]))

moon = CelestialBody('Moon', 7.34767309e22, 4.9048695e3, 1737.1, dist2parent=384400.0,
                     orbital_period=29 * 86400 + 12 * 3600.0 + 44 * 60.0 + 2.8)
