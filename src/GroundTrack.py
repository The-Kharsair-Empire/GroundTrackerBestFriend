import os
import matplotlib.pyplot as plt
from src.OrbitTool import eci2ecef

import numpy as np


def groundtracks(positions, times, titles, start_time=None, cities=None, show_plot=True, filename='groundtracks.png', dpi=300):

    if cities is None:
        cities = ['New York', 'Cairo', 'Mumbai', 'Moscow', 'Tokyo', 'Tianjin']
    plt.figure(figsize=(16, 8))
    # plt.style.use('dark_background')
    plot_coastlines()
    plot_cities(cities)
    assert len(positions) == len(times) == len(titles)

    for i in range(len(positions)):
        plot_groundtracks(positions[i], times[i], titles[i], start_time)
    plt.grid(linestyle='dotted')
    plt.xlim([-180, 180])
    plt.ylabel([-90, 90])
    plt.xlabel('Longitude (degrees $^\circ$)')
    plt.ylabel('Latitude (degrees $^\circ$)')
    plt.legend()

    if show_plot:
        plt.show()
    else:
        plt.savefig(filename, dpi=dpi)


def plot_groundtracks(rs: np.ndarray, ts: np.ndarray, name, start_time=None):

    spherical_coords = get_groundtracks(rs, ts, start_time)

    # plot first point bigger than the rest, and x is long so index 1 first
    plt.plot(spherical_coords[0, 1], spherical_coords[0, 0], 'ro', label=name)
    plt.plot(spherical_coords[1:, 1], spherical_coords[1:, 0], 'ro', markersize=1)

    pass


# TODO: this is obviously wrong, the eci to ecef part is not wrong, and ecef to latlong, the text book formula does
#  not work
def get_groundtracks(rs: np.ndarray, ts: np.ndarray, start_time=None):
    # if start_date is None, then the first time is at the start of the epoch time when theta gmt = 0
    # TODO: if not, read from spice file the earth theta gmt at given start time in gmt and start from that
    assert rs.shape[0] == ts.shape[0]
    day_in_sec = 86400
    # rs_ecef = np.zeros(rs.shape)
    latlongs = []
    for i in range(len(rs)):
        theta_gmt = ts[i] * np.pi * 2 / day_in_sec

        rs_ecef = eci2ecef(rs[i, :], theta_gmt)

        latlongs.append(rs_ecef2latlong(rs_ecef))
        # latlongs.append(rs_ecef2latlong(rs[i]))

    return np.array(latlongs)


def rs_ecef2latlong(r):
    r_mag = np.linalg.norm(r)
    # lat = np.arcsin(r[2] / r_mag)
    # long = np.arccos((r[0] / r_mag * np.cos(lat)))
    lat = np.arcsin(r[2] / r_mag)
    long = np.arctan2(r[1], r[0])

    # if r[1] / r_mag <= 0:
    #     long = 2 * np.pi - long
    r2d = 180.0 / np.pi

    return np.array([lat * r2d, long * r2d, r_mag])


def plot_coastlines():
    coastlines_coords_file = os.path.join(os.path.dirname(__file__), '..', 'file', 'csv', 'coastlines.csv')
    coastlines_coords = np.genfromtxt(coastlines_coords_file, delimiter=',')
    plt.plot(coastlines_coords[:, 0], coastlines_coords[:, 1], 'mo', markersize=0.3)
    return plt


def plot_cities(cities_to_plot=None):
    if cities_to_plot is None:
        cities_to_plot = []
    cities = get_world_city_coor()
    for n, city in enumerate(cities_to_plot):
        coords = cities[city]
        plt.plot([coords[1]], [coords[0]], 'ro', markersize=2)
        if n % 2 == 0:
            xytext = (0, 2)
        else:
            xytext = (0, -8)
        plt.annotate(city, [coords[1], coords[0]],
                     textcoords='offset points', xytext=xytext,
                     ha='center', color='y', fontsize='small')


def get_world_city_coor():
    data_file = os.path.join(os.path.dirname(__file__), '..', 'file', 'csv', 'world_cities.csv')
    with open(data_file, encoding="utf8") as f:
        lines = f.readlines()

    header = lines[0]
    cities = {}
    for line in lines[1:]:
        line = line.split(',')

        try:
            cities[line[1]] = (float(line[2]), float(line[3]))
        except Exception as e:
            print(e)

    return cities
