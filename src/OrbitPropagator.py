import numpy as np
from scipy.integrate import ode
from multiprocessing import Process, SimpleQueue
from src import CelestialData as celestialBody
from .TLE import tle2coes
from .OrbitTool import coes2rv, rv2coes, get_escape_velocity
from .Aerodynamics import calc_atm_density
from .SPICE import spice_load_kernels, spice_get_ephemeris_data
import spiceypy as spice


class OrbitPropagator:

    def __init__(self, r0, v0, timespan, timestep, body=celestialBody.earth, initial_mass=10.0):

        self.coes_rel = None
        self.current_step = 0
        self.main_body_spice_ephemeris = None
        self.spice_tspan = None
        self.n_bodies_ephemeris = None
        self.spice_start_time = None
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
        self.cache = {}

        self.perturbation = {
            'Thrust': False
        }
        if body == celestialBody.earth:
            self.perturbation['J2'] = False
            self.perturbation['Aero'] = False
            self.perturbation['N_bodies'] = False
            self.perturbation['SRP'] = False
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
                if arg == 'N_bodies':

                    spice_files = []
                    # list of dictionary of {name -> str: [data -> CelestialBody, spice_file -> str]}
                    self.pert_params['other_bodies'] = kwargs.get('other_bodies')
                    # self.pert_params['srp'] = kwargs.get('srp', False)  # Boolean
                    self.pert_params['frame'] = kwargs.get('frame')
                    # if self.pert_params['srp']:
                    #     self.pert_params['spice_file'] = kwargs.get('spice_file')  # parent body spice file
                    #     spice_files.append(self.pert_params['spice_file'])

                    for each_body in self.pert_params['other_bodies']:
                        spice_files.append(self.pert_params['other_bodies'][each_body][1])
                    spice_load_kernels(spice_files)
                    self.spice_start_time = spice.utc2et(kwargs.get('date_0', '2020-12-3'))
                    self.spice_tspan = np.linspace(self.spice_start_time, self.spice_start_time + self.timespan,
                                                   self.n_steps)
                    self.n_bodies_ephemeris = {
                        each_body: spice_get_ephemeris_data(each_body, self.spice_tspan,
                                                            self.pert_params['frame'], self.body.name)
                        for each_body in self.pert_params['other_bodies']
                    }
                    # if self.pert_params['srp']:
                    #     self.main_body_spice_ephemeris = spice_get_ephemeris_data(self.body.name, self.spice_tspan,
                    #                                                               self.pert_params['frame'], 'SUN')

                if arg == 'SRP':
                    # TODO: some unexpected results in plot, maybe the spice time isn't right?
                    spice_files = []
                    self.pert_params['mass'] = kwargs.get('mass', 10.0)
                    self.pert_params['CR'] = kwargs.get('CR')  # coefficient of reflection
                    self.pert_params['A_srp'] = kwargs.get('A_srp')  # area, for area to mass ratio
                    self.pert_params['G1'] = kwargs.get('G1')  # related to solar radiation intensity
                    self.pert_params['spice_file'] = kwargs.get('spice_file')  # parent body spice file
                    spice_files.append(self.pert_params['spice_file'])
                    spice_load_kernels(spice_files)

                    self.spice_start_time = spice.utc2et(kwargs.get('date_0', '2020-12-3'))
                    self.spice_tspan = np.linspace(self.spice_start_time, self.spice_start_time + self.timespan,
                                                   self.n_steps)
                    self.main_body_spice_ephemeris = spice_get_ephemeris_data(self.body.name, self.spice_tspan,
                                                                              kwargs.get('frame'), 'SUN')

    def calculate_all_coes(self, deg=True, ta_in_time=False, t=None):
        self.coes = np.zeros((self.n_steps, 6))
        for n in range(self.n_steps):
            self.coes[n, :] = rv2coes(self.rs[n, :], self.vs[n, :], self.body.mu, deg=deg, ta_in_time=ta_in_time, t=t)

        self.coes_rel = self.coes - self.coes[0, :]

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

        if 'J2' in self.perturbation and self.perturbation['J2']:
            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)

            a_j2 = (1.5 * self.body.j2 * self.body.mu * self.body.radius ** 2) / (norm_r ** 4) * np.array([tx, ty, tz])
            a += a_j2

        if 'Aero' in self.perturbation and self.perturbation['Aero']:
            # TODO: z should be geodetic altitude, this is geocentric, modify in the future
            z = norm_r - self.body.radius
            rho = calc_atm_density(z, self.body)
            v_rel = v - np.cross(self.body.angular_velocity, r)

            a_drag = -v_rel * np.linalg.norm(v_rel) * rho * self.pert_params['A'] * \
                     self.pert_params['Cd'] / (2 * m)

            a += a_drag

        if 'Thrust' in self.perturbation and self.perturbation['Thrust']:
            a += self.pert_params['direction'] * (v / np.linalg.norm(v)) * self.pert_params[
                'thrust'] / m / 1000.0  # km / s ** 2
            mdot = -self.pert_params['thrust'] / self.pert_params['isp'] / 9.81

        if 'N_bodies' in self.perturbation and self.perturbation['N_bodies']:
            for each_body in self.pert_params['other_bodies']:
                # vector from central body to nth perturbing body
                r_cb2body = self.n_bodies_ephemeris[each_body][self.current_step, :3]

                # vector from orbiting object to the nth perturbing body
                r_sat2body = r_cb2body - r

                nth_body_acc = self.pert_params['other_bodies'][each_body][0].mu * \
                               (r_sat2body / np.linalg.norm(r_sat2body) ** 3 - r_cb2body / np.linalg.norm(
                                   r_cb2body) ** 3)

                a += nth_body_acc

        if 'SRP' in self.perturbation and self.perturbation['SRP']:
            r_sun2sat = self.main_body_spice_ephemeris[self.current_step, :3] + r

            solar_rad_acc = (1 + self.pert_params['CR']) * self.pert_params['G1'] * self.pert_params['A_srp'] \
                            / m / np.linalg.norm(r_sun2sat) ** 3 * r_sun2sat

            a += solar_rad_acc

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

        self.current_step = 1

        solver = ode(self.ode_acc)
        solver.set_integrator(integrator)
        solver.set_initial_value(y0, 0)
        solver.set_f_params(self.body.mu)

        while solver.successful() and self.current_step < self.n_steps:
            solver.integrate(solver.t + self.dt)
            self.ts[self.current_step] = solver.t
            ys[self.current_step] = np.array(solver.y)
            self.current_step += 1
            if self.check_stop_condition(ys, self.current_step):
                break

        self.rs = ys[:self.current_step, :3]
        self.vs = ys[:self.current_step, 3:6]
        self.masses = ys[:self.current_step, 6]
        self.ts = self.ts[:self.current_step]
        self.n_steps = self.current_step

    def check_stop_condition(self, ys, current_step) -> bool:
        if np.linalg.norm(ys[current_step - 1, :3]) < self.body.radius:
            print("trajectory hit parent body surface")
            return True
        if 'min_alt' in self.stop_conditions.keys():
            min_alt = self.stop_conditions['min_alt'].get('min_alt')
            if np.linalg.norm(ys[current_step - 1, :3]) < min_alt + self.body.radius:
                print(f"altitude lower than minimum altitude {min_alt}")
                return True
        if 'max_alt' in self.stop_conditions.keys():
            max_alt = self.stop_conditions['max_alt'].get('max_alt')
            if np.linalg.norm(ys[current_step - 1, :3]) > max_alt + self.body.radius:
                print(f"altitude exceed msc altitude {max_alt}")
                return True
        if 'escape_velocity_reached' in self.stop_conditions.keys():
            current_r = ys[current_step - 1, :3]
            current_v = ys[current_step - 1, 3:6]
            escape_v = get_escape_velocity(current_r, self.body.mu)
            current_r_norm = np.linalg.norm(current_v)
            if current_r_norm > escape_v:
                print(f"escape velocity {escape_v} around {self.body.name} "
                      f"at height {current_r_norm} reached")
                return True
        if 'new_soi' in self.stop_conditions.keys():
            pass
        if 'min_mass' in self.stop_conditions.keys():
            pass
        if 'more' in self.stop_conditions.keys():
            pass
        return False

    def enable_stop_conditions(self, name, **kwargs):
        self.stop_conditions[name] = kwargs

    def get_time_stamps_1d(self):
        return self.ts.T[0]


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
