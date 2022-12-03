import numpy as np
from src import plot_n_orbit_3d
from src import CelestialData as cd
from src import OrbitPropagator
from src import coes2rv


def topic_1():
    body = cd.earth
    r_mag = body.radius + 500
    v_mag = np.sqrt(body.mu / r_mag)  # vis-viva equation for circular orbit

    r0 = [r_mag, 0, 0]
    v0 = [0, v_mag * 1.39, 0]  # 1.39 TLI ((7800 + 3119) / 7800), 1.41 Escape ((7800 + 3119 + 100) / 7800)

    timespan = 3600 * 24 * 6 * 3  # 3 day

    dt = 100.0  # timestep 100s

    propagator_1 = OrbitPropagator(r0, v0, timespan, dt, body=body)
    propagator_1.propagate_orbit('lsoda')

    body = cd.earth
    r_mag = body.radius + 500  # circular orbit at 500 km height
    v_mag = np.sqrt(body.mu / r_mag)  # vis-viva equation for circular orbit

    r0 = [r_mag, 0, 0]
    v0 = [0, v_mag, 0]  # right-handed frame, velocity is tangent to position at periapsis or apoapsis

    timespan = 100 * 60.0  # 90 mins, enough for one orbit in LEO

    dt = 100.0  # timestep 100s

    propagator_2 = OrbitPropagator(r0, v0, timespan, dt, body=body)
    propagator_2.propagate_orbit('lsoda')

    plot_n_orbit_3d([propagator_1.rs, propagator_2.rs], ['TLI', '500km-LEO'], body.radius)


def topic_2():
    body = cd.earth

    # ISS
    coes0 = (body.radius + 410.0, 0.0006189, 51.6393, 105.6372, 234.1955, 0.0)
    r0, v0 = coes2rv(*coes0, body.mu, True)

    # GEO sat
    coes1 = (body.radius + 35800, 0.0, 0.0, 75.0, 0.0, 0.0)
    r1, v1 = coes2rv(*coes1, body.mu, True)

    propagator_1 = OrbitPropagator(r0, v0, 60 * 100, 100.0, body)
    propagator_1.propagate_orbit('lsoda')

    propagator_2 = OrbitPropagator(r1, v1, 3600 * 24 * 6, 100.0, body)
    propagator_2.propagate_orbit('lsoda')

    plot_n_orbit_3d([propagator_1.rs, propagator_2.rs], ['ISS', 'GEO-sat'], body.radius)


if __name__ == '__main__':
    # topic_1()
    topic_2()

