from urllib.request import urlopen
import requests
import json
import pprint
from enum import Enum


# TODO: build this https api interface


class URLEncoding(Enum):
    SPACE = '%20'
    SEMICOLON = '%3B'


class CoordinateCenter(Enum):
    SOLAR_SYSTEM_BARYCENTRE = '@0'


class OutputUnit(Enum):
    KMS = 'KM-S'
    AUD = 'AU-D'
    KMD = 'KM-D'


class HORIZONSHttpClient:
    pass


def get_starman():  # just a test function of the horizon systems
    api_url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    params = {
        'format': 'text',
        'EPHEM_TYPE': 'VECTORS',
        'STEP_SIZE': '1d',
        'COMMAND': '-143205',
        'START_TIME': '2018-02-08',
        'STOP_TIME': '2018-02-09',
        'CENTER': '@0',
        'REF_PLANE': 'ECLIPTIC',
        'OUT_UNITS': 'KM-S',
        'REF_SYSTEM': 'ICRF',  # equivalent to J2000 in JPL ephemerides
    }
    response = requests.get(url=api_url, params=params)
    # page = urlopen(api_url)
    # html = page.read().decode("utf-8")
    # parsed = json.loads(html)
    # print(response.text)
    resp_lines = response.text.split('\n')

    start_from_line = 0
    for i in range(len(resp_lines)):
        # print(resp_lines[i])

        if resp_lines[i].startswith('$$SOE'):
            start_from_line = i + 1
            break

    data = []
    while not resp_lines[start_from_line].startswith('$$EOE'):
        line = resp_lines[start_from_line]
        if 'TDB' in line:
            data_line_1 = resp_lines[start_from_line + 1]
            data_line_2 = resp_lines[start_from_line + 2]
            data_line_3 = resp_lines[start_from_line + 3]
            print(resp_lines[start_from_line])
            print(data_line_1)
            print(data_line_2)
            print(data_line_3)
            X_start = data_line_1.find('X')
            Y_start = data_line_1.find('Y')
            Z_start = data_line_1.find('Z')
            VX_start = data_line_2.find('VX')
            VY_start = data_line_2.find('VY')
            VZ_start = data_line_2.find('VZ')

            X = data_line_1[X_start:Y_start]
            Y = data_line_1[Y_start:Z_start]
            Z = data_line_1[Z_start:]

            VX = data_line_2[VX_start:VY_start]
            VY = data_line_2[VY_start:VZ_start]
            VZ = data_line_2[VZ_start:]

            X = X[3:].strip()
            Y = Y[3:].strip()
            Z = Z[3:].strip()
            VX = VX[3:].strip()
            VY = VY[3:].strip()
            VZ = VZ[3:].strip()

            data.append({
                'x': float(X),
                'y': float(Y),
                'z': float(Z),
                'vx': float(VX),
                'vy': float(VY),
                'vz': float(VZ),
            })

        start_from_line += 1

    return data


if __name__ == '__main__':
    pprint.pprint(get_starman())
