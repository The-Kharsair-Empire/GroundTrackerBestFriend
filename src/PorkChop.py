# TODO: calculating porkchop plot value using parallel processing.
from .OrbitTool import coes2rv


def porkchop_plot_spice(planet0, planet1, center_body,
                        departure0, departure1, arrival0, arrival1,
                        search_grid_timestep_departure, search_grid_timestep_arrival,
                        reference_frame, title=None, config=None):
    _config = {}

    if title is None:
        title = f"{planet0} to {planet1} porkchop plot"


def porkchop_plot_generic(coes0, coes1, center_body,
                          departure_search_start,  # how many seconds since the time expressed in ta in coes0 should
                          # the departure search start
                          departure_window_width, departure_search_grid_timestep,
                          time_between_departure_search_start_to_arrival_search_start,  # first date on the x-axis to
                          # the first date on the y-axis
                          arrival_window_width, arrival_search_grid_timestep,
                          title=None):
    """
    time all in seconds
    """
    if title is None:
        title = "untitled porkchop plot"

    departure_search_count = departure_window_width // departure_search_grid_timestep
    arrival_search_count = arrival_window_width // arrival_search_grid_timestep

    # TODO: propagate the orbit of the planets from departure start to arrival end and acquire the whole ephemeris,
    # TODO:     the step size should be the search window size.

    # parallelization strategy: passing the range to find for each process,
    # i.e. process 1 find the first 10 search grid,
    # process 2 find the second 10,
    # so parameter is only the search range
    for i in range(departure_search_count):
        # r0, _ = coes2rv(coes0[0])
        for j in range(arrival_search_count):
            pass
