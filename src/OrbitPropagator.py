import numpy as np
from scipy.integrate import ode
from multiprocessing import Process, SimpleQueue
from src import CelestialData as celestialBody
from .TLE import tle2coes
from .Orbit import coes2rv


class OrbitPropagator:

    def __init__(self, r0, v0, timespan, timestep, body=celestialBody.earth):

        self.r0 = r0
        self.v0 = v0
        self.timespan = timespan
        self.dt = timestep
        self.body = body
        self.n_steps = int(np.ceil(timespan / timestep))
        self.vs = None
        self.rs = None
        self.perturbation = {
            'J2': False,
            'Aero': False,
            'Moon_Grav': False
        }

    def set_perturbation(self, *args):
        for arg in args:
            if arg in self.perturbation:
                self.perturbation[arg] = True

    def ode_two_body_acc(self, t, y, mu):
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])

        norm_r = np.linalg.norm(r)

        a = -mu * r / norm_r ** 3

        if self.perturbation['J2']:
            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)

            a_j2 = (1.5 * self.body.j2 * self.body.mu * self.body.radius ** 2) / (norm_r ** 4) * np.array([tx, ty, tz])
            a += a_j2

        ax, ay, az = a
        return [vx, vy, vz, ax, ay, az]

    def propagate_orbit(self, integrator='lsoda'):
        ys = np.zeros((self.n_steps, 6))
        ts = np.zeros((self.n_steps, 1))

        if isinstance(self.r0, np.ndarray):
            y0 = self.r0.tolist() + self.v0.tolist()
        else:
            y0 = self.r0 + self.v0

        ys[0] = np.array(y0)  # set initial state

        current_step = 1

        solver = ode(self.ode_two_body_acc)
        solver.set_integrator(integrator)
        solver.set_initial_value(y0, 0)
        solver.set_f_params(self.body.mu)

        while solver.successful() and current_step < self.n_steps:
            solver.integrate(solver.t + self.dt)
            ts[current_step] = solver.t
            ys[current_step] = np.array(solver.y)
            current_step += 1

        self.rs = ys[:, :3]
        self.vs = ys[:, 3:]


def task(r, v,
         positions: SimpleQueue, names: SimpleQueue, body, time, dt, solver):
    propagator = OrbitPropagator(r, v, time, dt, body)
    propagator.propagate_orbit(solver)
    positions.put(propagator.rs)
    # ts.put(tle.satellite_name)


def simulate_orbit_concurrently(list_of_parsed_tle, timespan, timestep, body=celestialBody.earth, solver='lsoda'):
    positions = SimpleQueue()
    names = SimpleQueue()
    processes = []
    for tle in list_of_parsed_tle:
        print(tle.satellite_name)
        print(f'inclination: {tle.inclination}')
        r, v = coes2rv(*(tle2coes(tle, body.mu)[:-1]), body.mu)
        p = Process(target=task, args=(r, v, positions, names, body, timespan, timestep, solver,))
        p.start()
        processes.append(p)

    rs = []
    titles = []

    for p in processes:
        p.join()
        # rs.append(positions.get())
        # titles
    # while not positions.empty():
    #     rs.append(positions.get())
    # titles.append(names.get())

    return rs, titles
