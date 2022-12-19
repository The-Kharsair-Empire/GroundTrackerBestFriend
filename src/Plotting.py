from typing import List

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

matplotlib.use('TkAgg')
plt.style.use('dark_background')

color_list: list[str] = ['r', 'y', 'b', 'm', 'g', 'c']
label_coloc_list_1 = ['b*', 'g*', 'r*', 'c*', 'y*', 'm*']
label_coloc_list_2 = ['c*', 'm*', 'y*', 'b*', 'r*', 'g*']


def plot_3d(r, body_rad, title='Orbit', show_plot=True, dpi=300):
    # 3D plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # plot central body
    r_plot = body_rad
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v)
    _y = r_plot * np.sin(_u) * np.sin(_v)
    _z = r_plot * np.cos(_v)
    ax.plot_wireframe(_x, _y, _z, cmap='Blues')

    # plot x, y, z vectors for reference
    l: float = r_plot * 2.0
    x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    u, v, w = [[l, 0, 0], [0, l, 0], [0, 0, l]]
    ax.quiver(x, y, z, u, v, w, color='k')

    # plot trajectory and starting point
    ax.plot(r[:, 0], r[:, 1], r[:, 2], 'r', label='Trajectory')
    ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'b*', label='Initial Position')

    # check for custom axes limits
    max_val = np.max(np.abs(r))

    # set Labels and title
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel('Z (km) ')

    # make x, y, z projected equal
    ax.set_aspect('equal')

    # ax.plot title()
    ax.set_title(title, color='y')
    plt.legend()
    # plt.legend(['Trajectory', 'Starting Position'])

    if show_plot:
        plt.show()
    else:
        plt.savefig('figures/{}.png'.format(title), dpi=dpi)


def plot_n_orbit_3d(rs, labels, body_rad, wire_frame=True, title='Orbits', highlight_endpoints=True, show_plot=True,
                    col=None, label_col_1=None, label_col_2=None, dpi=300):
    if label_col_1 is None:
        label_col_1 = label_coloc_list_1
    if label_col_2 is None:
        label_col_2 = label_coloc_list_2
    if col is None:
        col = color_list
    assert len(rs) == len(labels), "you should assign a label for each orbit"

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    r_plot = body_rad
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v)
    _y = r_plot * np.sin(_u) * np.sin(_v)
    _z = r_plot * np.cos(_v)
    if wire_frame:
        ax.plot_wireframe(_x, _y, _z, cmap='Blues')
    else:
        ax.plot_surface(_x, _y, _z, cmap='Blues')

    # plot x, y, z vectors for reference
    l: float = r_plot * 2.0
    x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    u, v, w = [[l, 0, 0], [0, l, 0], [0, 0, l]]
    ax.quiver(x, y, z, u, v, w, color='k')

    # plot trajectory and starting point
    max_val = 0
    for i in range(len(rs)):
        r = rs[i]
        label = labels[i]
        ax.plot(r[:, 0], r[:, 1], r[:, 2], col[i % len(col)], label=f'{label}\'s Trajectory')
        if highlight_endpoints:
            ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], label_col_1[i % len(label_col_1)],
                    label=f'{label}\'s Initial Position')
            ax.plot([r[-1, 0]], [r[-1, 1]], [r[-1, 2]], label_col_2[i % len(label_col_2)],
                    label=f'{label}\'s Final Position')
        this_max_val = np.max(np.abs(r))
        if this_max_val > max_val:
            max_val = this_max_val

    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel('Z (km) ')

    ax.set_aspect('equal')

    ax.set_title(title, color='y')
    plt.legend()

    if show_plot:
        plt.show()
    else:
        plt.savefig('figures/{}.png'.format(title), dpi=dpi)


def plot_coes_over_time(coes_array, ts, title='COEs vs Time', time_unit='hour', deg=True, show_plot=True,
                        figsize=(16, 8), dpi=300):
    if time_unit not in ['hour', 'day', 'second', 'minute']:
        raise ValueError("Supplied time_unit must be in one of['hour', 'day', 'second', 'minute']")
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=figsize)
    fig.suptitle(title, fontsize=20)

    # x axis
    if time_unit == 'hour':
        ts /= 3600.0
    elif time_unit == 'minute':
        ts /= 60.0
    elif time_unit == 'day':
        ts /= 3600.0 * 24

    if deg:
        angle_unit = 'deg'
    else:
        angle_unit = 'rad'

    axs[0, 0].plot(ts, coes_array[:, 5])
    axs[0, 0].set_title('True Anomaly ({}) vs Time ({})'.format(angle_unit, time_unit))
    axs[0, 0].grid(True)

    axs[0, 1].plot(ts, coes_array[:, 3])
    axs[0, 1].set_title('RAAN ({}) vs Time ({})'.format(angle_unit, time_unit))
    axs[0, 1].grid(True)

    axs[0, 2].plot(ts, coes_array[:, 4])
    axs[0, 2].set_title('Argument of Periapsis ({}) vs Time ({})'.format(angle_unit, time_unit))
    axs[0, 2].grid(True)

    axs[1, 0].plot(ts, coes_array[:, 0])
    axs[1, 0].set_title('Semi Major Axis (km) vs Time ({})'.format(time_unit))
    axs[1, 0].grid(True)

    axs[1, 1].plot(ts, coes_array[:, 1])
    axs[1, 1].set_title('Eccentricity vs Time ({})'.format(time_unit))
    axs[1, 1].grid(True)

    axs[1, 2].plot(ts, coes_array[:, 2])
    axs[1, 2].set_title('Inclination ({}) vs Time ({})'.format(angle_unit, time_unit))
    axs[1, 2].grid(True)

    if show_plot:
        plt.show()
    else:
        plt.savefig('figures/{}.png'.format(title), dpi=dpi)


def plot_one_parameter(ys, ts, ylabel, title='whatever', time_unit='hour', show_plot=True, figsize=(16, 8)
                            , dpi=300):
    if time_unit == 'hour':
        ts /= 3600.0
    elif time_unit == 'minute':
        ts /= 60.0
    elif time_unit == 'day':
        ts /= 3600.0 * 24

    plt.figure(figsize=figsize)
    plt.plot(ts, ys, 'w')
    plt.grid(True)
    plt.xlabel(time_unit)
    plt.ylabel(ylabel)
    plt.title(title)

    if show_plot:
        plt.show()
    else:
        plt.savefig('figures/{}.png'.format(title), dpi=dpi)


def plot_altitude_over_time(rs, ts, body, title="altitude over time", time_unit='hour', show_plot=True, figsize=(16, 8)
                            , dpi=300):
    altitudes = np.linalg.norm(rs, axis=1) - body.radius
    if time_unit == 'hour':
        ts /= 3600.0
    elif time_unit == 'minute':
        ts /= 60.0
    elif time_unit == 'day':
        ts /= 3600.0 * 24

    plt.figure(figsize=figsize)
    plt.plot(ts, altitudes, 'w')
    plt.grid(True)
    plt.xlabel(time_unit)
    plt.ylabel('altitude (km)')
    plt.title(title)

    if show_plot:
        plt.show()
    else:
        plt.savefig('figures/{}.png'.format(title), dpi=dpi)


def plot_apoapsis_n_periapsis_over_time(aps, pes, ts, body, title="Ap, Pe over time", time_unit='hour', show_plot=True
                                        , figsize=(16, 8)
                                        , dpi=300):
    if time_unit == 'hour':
        ts /= 3600.0
    elif time_unit == 'minute':
        ts /= 60.0
    elif time_unit == 'day':
        ts /= 3600.0 * 24

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=figsize)
    fig.suptitle(title, fontsize=20)

    max_scale = np.max(aps)
    min_scale = np.min(pes)

    axs[0].plot(ts, aps)
    axs[0].set_title('Apoapsis (km) vs Time ({})'.format(time_unit))
    axs[0].grid(True)
    axs[0].set_ylim([min_scale, max_scale])

    axs[1].plot(ts, pes)
    axs[1].set_title('Periapsis (km) vs Time ({})'.format(time_unit))
    axs[1].grid(True)
    axs[1].set_ylim([min_scale, max_scale])

    if show_plot:
        plt.show()
    else:
        plt.savefig('figures/{}.png'.format(title), dpi=dpi)
