from dataclasses import dataclass


@dataclass
class CelestialBody:
    name: str
    mass: float  # kg
    mu: float  # km^3 / s^2
    radius: float  # km
    j2: float


sun = CelestialBody('Sun', 1.989e30, 1.32712e11, 695700.0, 0.0)
earth = CelestialBody('Earth', 5.972e24, 3.986e5, 6378.0, -0.001082635854)
