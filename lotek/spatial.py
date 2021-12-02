"""
Lotek processing helper functions for spatial analyses
"""


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
    dist_tmp = np.linalg.norm(start_coords-end_coords, axis=1)
    return dist_tmp

