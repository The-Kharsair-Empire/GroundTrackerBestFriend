import numpy as np


def calc_atm_density(z, body):
    rhos, zs = None, None
    if z > 1000.0:  # assume no drag 1000km above
        rhos, zs = (0.0, 0.0), (0.0, 0.0)
    else:
        assert len(body.rhos) == len(body.zs)
        # assert z > body.zs[0], f"altitude {z} below base reference density {body.zs[0]}"
        for n in range(len(body.rhos) - 1):
            if body.zs[n] < z < body.zs[n + 1]:
                rhos = (body.rhos[n], body.rhos[n + 1])
                zs = (body.zs[n], body.zs[n + 1])
                break

    if rhos is None:
        # raise Exception("Cannot get reference density")
        rhos, zs = (0.0, 0.0), (0.0, 0.0)
    if rhos[0] == 0.0:
        return 0.0

    scale_height = -(zs[1] - zs[0]) / np.log(rhos[1] / rhos[0])

    return rhos[0] * np.exp(-(z - zs[0]) / scale_height)
