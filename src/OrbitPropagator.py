import numpy as np
from scipy.integrate import ode

from src import CelestialData as celestialBody


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

    @staticmethod
    def ode_two_body_acc(t, y, mu):
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])

        norm_r = np.linalg.norm(r)

        ax, ay, az = -mu * r / norm_r ** 3

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
