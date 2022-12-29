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


def make_steps(x, start_field_x='UTM_X_fnl_lag1', start_field_y='UTM_Y_fnl_lag1',
               end_field_x='UTM_X_fnl', end_field_y='UTM_Y_fnl'):
    from shapely.geometry import Point, LineString
    """
    Create line-strings (i.e., polylines) between each GPS start/end
    ***NOTE: intended to be used with pandas.DataFrame.apply()***
    
    Parameters
    -----------
    start_field_x :      (string) name of the field with starting x coordinates
    start_field_y :      (string) name of the field with starting y coordinates
    end_field_x :      (string) name of the field with ending x coordinates
    end_field_y :      (string) name of the field with ending y coordinates
    
    Returns
    -----------
    LineString() object from shapely module
    """
    line_step = LineString([Point(x[start_field_x], x[start_field_y]),
                 Point(x[end_field_x], x[end_field_y])])
    return line_step


def graze_intensity(x):
    """
    Calculate average grazing intensity for a line segment connecting two GPS locations
    based on the duration ('Fix_Duration') of time and the distance in meters ('steplength') 
    between the two fixes
    ***NOTE: intended to be used with geopandas.GeoDataFrame.apply()***
    
    Returns
    ---------
    A float value representing the seconds spent grazing per meter along the line
    """
    
    if x['steplength'] == 0 or x['cell_ct'] == 0:
        return 60.0*x['Fix_Duration']
    else:
        gi = (60.0*x['Fix_Duration'])/x['cell_ct']
    if gi > 300:
        return 60.0*x['Fix_Duration']
    else:
        return gi