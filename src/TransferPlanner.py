from .OrbitTool import hohmann_transfer_scalar_dv
from .Plotting import plot_n_orbit_3d


def transfer_planner(coes0, coes1, body):
    # TODO:
    """
    assumption, co-planar, impulsive burn at pe or ap only.
    all orbit are elliptical.
    spacecraft in prograde direction

    """

    a_init, e_init, i_init, raan_init, aop_init, ta_init = coes0
    a_final, e_final, i_final, raan_final, aop_final, ta_final = coes1

    ap_init = a_init * (1 - e_init)
    ap_final = a_final * (1 - e_final)

    pe_init = a_init * (1 + e_init)
    pe_final = a_final * (1 + e_final)



