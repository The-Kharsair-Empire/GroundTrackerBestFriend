import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

matplotlib.use('TkAgg')
plt.style.use('dark_background')


def plot_3d(r, body_rad, title='Orbit', show_plot=True):
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
        plt.savefig('figures/{}.png'.format(title), dpi=300)


def plot_n_orbit_3d(rs, labels, body_rad, col=['r', 'y', 'b', 'm'], label_col=['b*', 'g*', 'r*', 'y*'], title='Orbits', show_plot=True):
    assert len(rs) == len(labels), "you should assign a label for each orbit"

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

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
    max_val = 0
    for i in range(len(rs)):
        r = rs[i]
        label = labels[i]
        ax.plot(r[:, 0], r[:, 1], r[:, 2], col[i], label=f'{label}\'s Trajectory')
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], label_col[i], label=f'{label}\'s Initial Position')
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
        plt.savefig('figures/{}.png'.format(title), dpi=300)