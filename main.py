import datetime

import numpy as np
from src import plot_n_orbit_3d \
    , plot_coes_over_time, plot_altitude_over_time \
    , plot_apoapsis_n_periapsis_over_time, plot_one_parameter
from src import CelestialData as cd
from src import OrbitPropagator, simulate_orbit_concurrently
from src import coes2rv
from src import tle2coes, parse_raw_tle
from src import spice_load_kernels, spice_get_spk_objects, tc2array, spice_get_ephemeris_data


def task_1():  # Two body equation of motion
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


def task_2():  # Classical Orbital Elements
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


def task_3():  # TLE
    tle_raw = \
        """
ISS (ZARYA)
1 25544U 98067A   22339.03402137  .00010856  00000+0  19968-3 0  9994
2 25544  51.6438 210.2304 0004354 127.5503 200.9328 15.49745241371670
CSS (TIANHE)
1 48274U 21035A   22338.24722390  .00020998  00000+0  26284-3 0  9995
2 48274  41.4771 263.3901 0005969  11.7317 116.1668 15.59272697 91354
CZ-6A DEB               
1 54265U 22151G   22338.72774229  .00000471  00000+0  63325-3 0  9994
2 54265  98.9964 341.9167 0384511  18.1840 343.2615 13.37528759  3012
ORION                   
1 54257U 22156A   22320.38000000 -.00015874  00000+0  00000+0 0  9993
2 54257  30.5280  10.5600 9646414  20.9260   1.0490  0.10079459    30
"""
    body = cd.earth
    tle_parsed = parse_raw_tle(tle_raw)
    rs = []
    titles = []

    start_linear = datetime.datetime.now()
    for tle in tle_parsed:
        print(tle.satellite_name)
        print(f'inclination: {tle.inclination}')
        coes = tle2coes(tle, body.mu)
        r, v = coes2rv(*coes[:-1], body.mu)
        propagator = OrbitPropagator(r, v, 40 * 3600 * 24, 100.0, body)
        propagator.propagate_orbit('lsoda')
        rs.append(propagator.rs)
        titles.append(tle.satellite_name)

    end_linear = datetime.datetime.now()
    print(f'time taken, linear {(end_linear - start_linear).total_seconds()}')

    # start_mp = datetime.datetime.now()
    # # motherfuck this doesn't work
    # rs, titles = simulate_orbit_concurrently(tle_parsed, 40 * 3600 * 24, 100.0, body, 'lsoda')
    #
    # end_mp = datetime.datetime.now()
    # print(f'time taken, multi-processing {(end_mp - start_mp).total_seconds()}')
    plot_n_orbit_3d(rs, titles, body.radius)


def task_4():  # CZ-6 breakup debris
    file = open('file/tle/CZ6A_DEB.txt', 'r')
    tle_raw = file.read()
    file.close()
    body = cd.earth
    tle_parsed = parse_raw_tle(tle_raw)
    rs = []
    titles = []

    for tle in tle_parsed:
        print(tle.satellite_name)
        print(f'inclination: {tle.inclination}')
        coes = tle2coes(tle, body.mu)
        r, v = coes2rv(*coes[:-1], body.mu)
        propagator = OrbitPropagator(r, v, 60 * 170, 100.0, body)
        propagator.propagate_orbit('lsoda')
        rs.append(propagator.rs)
        titles.append(tle.satellite_name)

    plot_n_orbit_3d(rs, titles, body.radius)


def task_5():  # J2 perturbation
    tle_raw = \
        """
ISS (ZARYA)
1 25544U 98067A   22339.03402137  .00010856  00000+0  19968-3 0  9994
2 25544  51.6438 210.2304 0004354 127.5503 200.9328 15.49745241371670
"""
    body = cd.earth
    tle_parsed = parse_raw_tle(tle_raw)
    rs = []
    titles = []

    for tle in tle_parsed:
        coes = tle2coes(tle, body.mu)
        r, v = coes2rv(*coes[:-1], body.mu)
        propagator = OrbitPropagator(r, v, 3600 * 24 * 3, 100.0, body)
        propagator.enable_perturbation('J2')
        propagator.propagate_orbit('lsoda')
        rs.append(propagator.rs)
        titles.append(tle.satellite_name)

    plot_n_orbit_3d(rs, titles, body.radius, True, 'J2 perturbation')


def task_6():  # coes from r, v state vector
    tle_raw = \
        """
ISS (ZARYA)
1 25544U 98067A   22339.03402137  .00010856  00000+0  19968-3 0  9994
2 25544  51.6438 210.2304 0004354 127.5503 200.9328 15.49745241371670
"""
    body = cd.earth
    tle_parsed = parse_raw_tle(tle_raw)
    coes = tle2coes(tle_parsed[0], body.mu)
    r, v = coes2rv(*coes[:-1], body.mu)
    propagator = OrbitPropagator(r, v, 3600 * 24 * 30, 100.0, body)
    propagator.enable_perturbation('J2')
    propagator.propagate_orbit('lsoda')
    propagator.calculate_all_coes(deg=True)

    plot_coes_over_time(propagator.coes, propagator.ts, time_unit='day')


def task_7():  # sun synchronous orbit

    body = cd.earth
    r, v = coes2rv(body.radius + 600, 0.01, 63.435, 0.0, 0.0, 50.0, body.mu, deg=True)
    propagator = OrbitPropagator(r, v, 3600 * 48, 100.0, body)
    propagator.enable_perturbation('J2')
    propagator.propagate_orbit('lsoda')
    propagator.calculate_all_coes(deg=True)

    plot_coes_over_time(propagator.coes, propagator.ts, time_unit='hour')


def task_8():  # air-drag

    body = cd.earth
    pe = 215 + body.radius
    ap = 300 + body.radius
    a = (pe + ap) / 2.0
    e = (ap - pe) / (ap + pe)
    raan = 340.0
    i = 62.5
    aop = 58.0
    ta = 350.0

    r, v = coes2rv(a, e, i, raan, aop, ta, body.mu, deg=True)
    propagator = OrbitPropagator(r, v, 3600 * 24 * 5, 100.0, body)
    # propagator.enable_perturbation('J2')
    propagator.enable_perturbation('Aero', Cd=2.2, A=(1e-3 ** 2 / 4.0), mass=10.0)
    propagator.propagate_orbit('lsoda')

    plot_n_orbit_3d([propagator.rs], ['Spiral Re-entry'], body.radius)

    plot_altitude_over_time(propagator.rs, propagator.ts, body)

    propagator.calculate_ap_pe()

    plot_apoapsis_n_periapsis_over_time(propagator.apoapsis, propagator.periapsis, propagator.ts, body)

    plot_coes_over_time(propagator.coes, propagator.ts, time_unit='hour')


def task_9():  # thrust trajectory

    body = cd.earth
    a = body.radius + 20000
    e = 0.6
    raan = 200.0
    i = 10.0
    aop = 30.0
    ta = 0.0

    r, v = coes2rv(a, e, i, raan, aop, ta, body.mu, deg=True)
    propagator = OrbitPropagator(r, v, 3600 * 24, 100.0, body)
    # propagator.enable_perturbation('J2')

    # thrust in newton, mass in kg
    propagator.enable_perturbation('Thrust', isp=4300.0, direction=1.0, thrust=0.127, mass=10.0)

    propagator.enable_perturbation('Aero', Cd=2.2, A=(1e-3 ** 2 / 4.0), mass=10.0)
    propagator.enable_perturbation('J2')
    # propagator.enable_stop_conditions('min_alt', min_alt=500)
    propagator.propagate_orbit('lsoda')

    plot_altitude_over_time(propagator.rs, propagator.ts, body)

    propagator.calculate_ap_pe()

    plot_apoapsis_n_periapsis_over_time(propagator.apoapsis, propagator.periapsis, propagator.ts, body)

    plot_coes_over_time(propagator.coes, propagator.ts, time_unit='hour')

    plot_n_orbit_3d([propagator.rs], ['low thrust trajectory'], body.radius, True, 'Orbit')
    # plot_one_parameter(propagator.masses, propagator.ts, 'mass', 'spacecraft mass change over time')


def task_10():  # visualize solar system from spice data

    frame = 'ECLIPJ2000'
    # frame = 'J2000' # Earth rotation axis as reference
    body = cd.sun
    observer = 'SUN'
    spice_load_kernels([
        'latest_leapseconds.tls.pc',
        'de440s.bsp'
    ])
    ids, names, tcs_sec, tcs_cal = spice_get_spk_objects('de440s.bsp', True)
    # ids, names, tcs_sec, tcs_cal = spice_get_spk_objects('de432s.bsp', True)
    times = tc2array(tcs_sec[0], 100000)
    names = list(filter(lambda x: 'BARYCENTER' in x, names))
    rs = []

    for name in names:
        rs.append(spice_get_ephemeris_data(name, times, frame, observer))

    plot_n_orbit_3d(rs, names, body.radius, False, 'solar system orbit', False)


def task_11():
    # TODO: propagate orbit for Pluto and Neptune in de432.bsp
    pass


def task_12():
    # TODO: n-body propagation

    body = cd.earth
    timespan = 3600 * 24 * 100.0
    # file = open('file/tle/ISS.txt', 'r')
    # tle_raw = file.read()
    # file.close()
    # iss_tle = parse_raw_tle(tle_raw)[0]
    # iss_coes = tle2coes(iss_tle, body.mu)
    # r, v = coes2rv(*iss_coes[:-1], body.mu)
    rs = []
    titles = ['GEO Sat', 'Moon']
    #
    frame = 'ECLIPJ2000'
    date_0 = '2020-02-23'
    # propagator = OrbitPropagator(r, v, timespan, 400.0, body)
    # propagator.enable_perturbation('N_bodies',
    #                                srp=True, frame=frame, spice_file='de440s.bsp',
    #                                date_0='2020-02-23',
    #                                other_bodies={
    #                                    'MOON': (cd.moon, 'de440s.bsp')
    #                                })
    # propagator.propagate_orbit('dopri5')
    # propagator.calculate_all_coes(deg=True)
    #
    # plot_coes_over_time(propagator.coes, propagator.ts, time_unit='hour')
    # rs.append(propagator.rs)

    a = 42164.0
    e = 0.001
    raan = 100.0
    i = 0.1
    aop = 200.0
    ta = 0.0

    r, v = coes2rv(a, e, i, raan, aop, ta, body.mu, deg=True)

    propagator = OrbitPropagator(r, v, timespan, 5000.0, body)
    propagator.enable_perturbation('N_bodies',
                                   frame=frame, spice_file='de440s.bsp',
                                   date_0=date_0,
                                   other_bodies={
                                       'MOON': (cd.moon, 'de440s.bsp')
                                   })
    propagator.propagate_orbit('dopri5')
    propagator.calculate_all_coes(deg=True)

    plot_coes_over_time(propagator.coes, propagator.ts, time_unit='day')
    rs.append(propagator.rs)
    rs.append(propagator.n_bodies_ephemeris['MOON'][:, :3])

    plot_n_orbit_3d(rs, titles, body.radius, True, 'comparing Moon perturbation to objects orbiting Earth', True)


if __name__ == '__main__':
    # task_1()
    # task_2()
    # task_3()
    # task_4()
    # task_5()
    # task_6()
    # task_7()
    # task_9()
    # task_10()
    task_12()
