import numpy as np


def calc_ta(a_list, b_list, c_list):
    """
    Calculate turning angle between fixes stored in three lists

    Paramters
    ---------
    a_list : list of previous fix coordinates (e.g., current fix lagged -1)
    b_list : list of current fix coordinates
    c_list : list of next fix coordinates (e.g., current fix lagged +1)

    Returns
    -------
    list of turning angles in degrees as departure from straight line heading (possible values: 0 - 180)
    """
    out_series = []
    for a, b, c in zip(a_list, b_list, c_list):
        ba = a - b
        bc = c - b

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        out_series.append(abs(180 - np.degrees(angle)))
    return out_series


def calc_dist(start_coords, end_coords):
    """
    Calculate distance between two lists of UTM coordinates

    Parameters
    ----------
    start_coords :     (list) starting UTM coordinates
    end_coords :       (list) ending UTM coordinates

    Returns
    -------
    list of distances between each coordinate pair in UTM units
    """
    dist_tmp = np.linalg.norm(start_coords - end_coords, axis=1)
    return dist_tmp