import numpy as np
from scipy.integrate import ode
from multiprocessing import Process, SimpleQueue
from src import CelestialData as celestialBody
from .TLE import tle2coes
from .OrbitTools import coes2rv, rv2coes
from .Aerodynamics import calc_atm_density


class OrbitPropagator:

    def __init__(self, r0, v0, timespan, timestep, body=celestialBody.earth, initial_mass=10.0):

        self.masses = None
        self.apoapsis = None
        self.periapsis = None
        self.r0 = r0
        self.v0 = v0
        self.timespan = timespan
        self.dt = timestep
        self.body = body
        self.n_steps = int(np.ceil(timespan / timestep))
        self.vs = None
        self.rs = None
        self.coes = None
        self.ts = None
        self.perturbation = {
            'Thrust': False
        }
        if body == celestialBody.earth:
            self.perturbation['J2'] = False
            self.perturbation['Aero'] = False
            self.perturbation['Moon_Grav'] = False
        self.pert_params = {
            'mass': initial_mass  # kg, spacecraft initial mass
        }

        self.stop_conditions = {}

    def enable_perturbation(self, *args, **kwargs):
        for arg in args:
            if arg in self.perturbation:
                self.perturbation[arg] = True
                if arg == 'Aero':
                    self.pert_params['Cd'] = kwargs.get('Cd')
                    self.pert_params['A'] = kwargs.get('A')
                    self.pert_params['mass'] = kwargs.get('mass')  # mass of the spacecraft, in kg
                if arg == 'Thrust':
                    self.pert_params['thrust'] = kwargs.get('thrust')
                    self.pert_params['isp'] = kwargs.get('isp')
                    self.pert_params['direction'] = kwargs.get('direction')
                    self.pert_params['mass'] = kwargs.get('mass')

    def calculate_all_coes(self, deg=True, ta_in_time=False, t=None):
        self.coes = np.zeros((self.n_steps, 6))
        for n in range(self.n_steps):
            self.coes[n, :] = rv2coes(self.rs[n, :], self.vs[n, :], self.body.mu, deg=deg, ta_in_time=ta_in_time, t=t)

    def calculate_ap_pe(self):
        if self.coes is None:
            self.calculate_all_coes()
        self.apoapsis = self.coes[:, 0] * (1 + self.coes[:, 1])
        self.periapsis = self.coes[:, 0] * (1 - self.coes[:, 1])

    def ode_acc(self, t, y, mu):
        rx, ry, rz, vx, vy, vz, m = y
        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])

        norm_r = np.linalg.norm(r)

        a = -mu * r / norm_r ** 3

        mdot = 0

        if self.perturbation['J2']:
            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)

            a_j2 = (1.5 * self.body.j2 * self.body.mu * self.body.radius ** 2) / (norm_r ** 4) * np.array([tx, ty, tz])
            a += a_j2

        if self.perturbation['Aero']:
            # TODO: z should be geodetic altitude, this is geocentric, modify in the future
            z = norm_r - self.body.radius
            rho = calc_atm_density(z, self.body)
            v_rel = v - np.cross(self.body.angular_velocity, r)

            a_drag = -v_rel * np.linalg.norm(v_rel) * rho * self.pert_params['A'] * \
                     self.pert_params['Cd'] / (2 * m)

            a += a_drag

        if self.perturbation['Thrust']:
            a += self.pert_params['direction'] * (v / np.linalg.norm(v)) * self.pert_params['thrust'] / m / 1000.0  # km / s ** 2
            mdot = -self.pert_params['thrust'] / self.pert_params['isp'] / 9.81
        ax, ay, az = a
        return [vx, vy, vz, ax, ay, az, mdot]

    def propagate_orbit(self, integrator='lsoda'):
        ys = np.zeros((self.n_steps, 7))
        self.ts = np.zeros((self.n_steps, 1))

        if isinstance(self.r0, np.ndarray):
            y0 = self.r0.tolist() + self.v0.tolist()
        else:
            y0 = self.r0 + self.v0

        y0 += [self.pert_params['mass']]

        ys[0] = np.array(y0)  # set initial state

        current_step = 1

        solver = ode(self.ode_acc)
        solver.set_integrator(integrator)
        solver.set_initial_value(y0, 0)
        solver.set_f_params(self.body.mu)

        while solver.successful() and current_step < self.n_steps:
            solver.integrate(solver.t + self.dt)
            self.ts[current_step] = solver.t
            ys[current_step] = np.array(solver.y)
            current_step += 1
            if self.check_stop_condition(ys, current_step):
                break

        self.rs = ys[:current_step, :3]
        self.vs = ys[:current_step, 3:6]
        self.masses = ys[:current_step, 6]
        self.ts = self.ts[:current_step]
        self.n_steps = current_step

    def check_stop_condition(self, ys, current_step) -> bool:
        if np.linalg.norm(ys[current_step-1, :3]) < self.body.radius:
            print("trajectory hit parent body surface")
            return True
        if 'min_alt' in self.stop_conditions.keys():
            min_alt = self.stop_conditions['min_alt'].get('min_alt')
            if np.linalg.norm(ys[current_step-1, :3]) < min_alt + self.body.radius:
                print(f"altitude lower than minimum altitude {min_alt}")
                return True
        if 'max_alt' in self.stop_conditions.keys():
            max_alt = self.stop_conditions['max_alt'].get('max_alt')
            if np.linalg.norm(ys[current_step-1, :3]) > max_alt + self.body.radius:
                print(f"altitude exceed msc altitude {max_alt}")
        if 'new_soi' in self.stop_conditions.keys():
            pass
        if 'min_mass' in self.stop_conditions.keys():
            pass
        if 'more' in self.stop_conditions.keys():
            pass
        return False

    def enable_stop_conditions(self, name, **kwargs):
        self.stop_conditions[name] = kwargs


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
