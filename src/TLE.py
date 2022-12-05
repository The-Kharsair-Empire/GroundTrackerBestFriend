import datetime
import numpy as np
import string
from dataclasses import dataclass

from .Orbit import M2E, E2ta
from .CelestialData import earth


@dataclass
class TLE:
    """
    TODO: not done, some of them in raw string, like second derivatives,
    TODO: need to convert them to proper float value
    """
    satellite_name: str
    satellite_number: int
    classification: str
    international_designator: str
    epoch_year: str
    epoch: str
    first_derivative: float
    second_derivative: str
    bstar: str
    ephemeris_type: int
    element_number: int
    inclination: float
    raan: float
    eccentricity: float
    argument_of_periapsis: float
    mean_anomaly: float
    mean_motion: float
    revolution_number: int


def checksum(line):
    # stolen from https://space.stackexchange.com/questions/5358/what-does-check-sum-tle-mean
    L = line.strip()
    cksum = 0
    for i in range(68):
        c = L[i]
        if c == ' ' or c == '.' or c == '+' or c in string.ascii_letters:
            continue
        elif c == '-':
            cksum = cksum + 1
        else:
            cksum = cksum + int(c)

    cksum %= 10

    return cksum


def parse_raw_tle(tle_text: str, perform_check=True):
    tle_lines = tle_text.split('\n')
    tle_lines = list(
        filter(
            lambda x: x != '',
            map(
                lambda x: x.strip('\t\n '),
                tle_lines)))

    TLEs = []
    for index in range(0, len(tle_lines), 3):
        name = tle_lines[index]
        line1 = tle_lines[index + 1]
        line2 = tle_lines[index + 2]

        # check sum
        if perform_check:
            check_sum_line_1 = int(line1[-1])
            check_sum_line_2 = int(line2[-1])
            if checksum(line1) != check_sum_line_1:
                raise Exception("line 1 checksum doesn't match!!")
            if checksum(line2) != check_sum_line_2:
                raise Exception("line 2 checksum doesn't match!!")

        tle = TLE(
            name, int(line1[2:7].strip()), line1[7], line1[9:17].strip(),
            line1[18:20].strip(), line1[20:32].strip(),
            float(line1[33:43].strip()), line1[44:52].strip(),
            line1[53:61].strip(), int(line1[62].strip()), int(line1[64:68].strip()),
            float(line2[8:16].strip()), float(line2[17:25].strip()), float("0." + line2[26:33].strip()),
            float(line2[34:42].strip()), float(line2[43:51].strip()), float(line2[52:63].strip()),
            int(line2[63:68].strip())
        )

        TLEs.append(tle)

    return TLEs


def tle2coes(tle: TLE, mu=earth.mu):

    yy, mm, dd, hh = calc_epoch(tle.epoch_year, tle.epoch)

    d2r = np.pi / 180.0
    i = tle.inclination * d2r
    raan = tle.raan * d2r
    aop = tle.argument_of_periapsis * d2r

    M = tle.mean_anomaly * d2r

    T = 1 / tle.mean_motion * 24 * 3600  # period in seconds

    a = (mu ** (1. / 3)) / ((2 * tle.mean_motion * np.pi) / 86400) ** (2. / 3)
    a_2 = (T ** 2 * mu / 4.0 / np.pi ** 2) ** (1 / 3.0)

    if not np.isclose(a, a_2):
        raise Exception("a doesn't match, check plugged in value or algorithm")

    E = M2E(M, tle.eccentricity)

    ta = E2ta(E, tle.eccentricity)

    return a, tle.eccentricity, i, raan, aop, ta, (yy, mm, dd, hh)


def calc_epoch(epoch_year: str, epoch: str):
    year = int('20' + epoch_year)

    day_of_the_year, day_fraction = epoch.split('.')
    day_of_the_year = int(day_of_the_year) - 1

    hour = float('0' + day_fraction) * 24

    date = datetime.date(year, 1, 1) + datetime.timedelta(day_of_the_year)

    month = float(date.month)
    day = float(date.day)

    return year, month, day, hour


