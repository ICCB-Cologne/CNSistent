import numpy as np


def count_below_lim(vals, min_val=0, max_val=1, steps=1000):
    """
    Counts how many samples are below each cutoff value, linearly spaced.

    Parameters
    ----------
    vals : array-like
        Array of values to count.
    min_val : float, optional
        Minimum value for the cutoff range. Default is 0.
    max_val : float, optional
        Maximum value for the cutoff range. Default is 1.
    steps : int, optional
        Number of steps between min_val and max_val. Default is 1000.

    Returns
    -------
    cutoffs : numpy.ndarray
        Array of cutoff values.
    counts : numpy.ndarray
        Array of counts of samples below each cutoff value.
    """
    cutoffs = np.linspace(min_val, max_val, steps + 1)
    #  Finds the indices where elements should be inserted to maintain order, 
    #  effectively counting the number of elements less than or equal to each cutoff.
    counts = np.searchsorted(np.sort(vals), cutoffs, side='right') / len(vals)	
    return cutoffs, counts


def count_cum_val(vals, min_val=0, max_val=1):
    """
    Counts the number of samples below each present value.

    Parameters
    ----------
    vals : array-like
        Array of values to count.
    min_val : float, optional
        Minimum value to consider. Default is 0.
    max_val : float, optional
        Maximum value to consider. Default is 1.

    Returns
    -------
    unique_vals : numpy.ndarray
        Array of unique values within the specified range.
    cumulative_count : numpy.ndarray
        Cumulative count of samples below each unique value.
    """
    vals = np.array(vals)
    vals = vals[(vals >= min_val) & (vals <= max_val)]
    unique_vals, counts = np.unique(vals.astype(np.float32), return_counts=True)
    cumulative_count = np.cumsum(counts) / len(vals)
    return unique_vals, cumulative_count


def calculate_signed_angle(s1, s2):
    """
    Calculates the signed angle between two slopes.
    Concave angles (rising curve) are positive, convex angles (falling curve) are negative.

    Parameters
    ----------
    s1 : float
        Slope of the first line.
    s2 : float
        Slope of the second line.

    Returns
    -------
    angle_degrees : float
        Signed angle in degrees.
    """
    angle_radians = np.arctan((s2 - s1) / (1 + s1 * s2))
    
    # Convert the angle to degrees
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees


def _get_direction(y):
    """
    Determines the direction of a monotonic curve.

    Parameters
    ----------
    y : array-like
        Array of y-values.

    Returns
    -------
    int
        1 if the curve is non-decreasing, -1 if the curve is non-increasing.

    Raises
    ------
    ValueError
        If the curve is non-monotonic.
    """
    if  np.all(np.diff(y) >= 0):
        return 1
    elif np.all(np.diff(y) < 0):
        return -1
    else:
        raise ValueError('Trying to find a knee in a curve that is non-monotonic')


def find_knee(x, y, knee=True):
    """
    Finds a knee or elbow in the curve using convex/concave curves.

    Parameters
    ----------
    x : array-like
        Array of x-values.
    y : array-like
        Array of y-values.
    knee : bool, optional
        If True, finds the knee. If False, finds the elbow. Default is True.

    Returns
    -------
    int
        Index of the knee or elbow.
    float
        Value of the knee or elbow.
    """
    if len(x) < 2:
        return -1, np.nan
    y_range = y[-1] - y[0]
    x_range = x[-1] - x[0]
    if y_range == 0:
        return -1, np.nan
    y = (np.array(y) - y[0]) / y_range 
    x = (np.array(x) - x[0]) / x_range
    orientation = _get_direction(y)
    orientation *= (-1 if knee else 1)

    # calculate the difference between slopes on the left and right side of each point 
    angles = []
    slopes = []    
    for i in range(1, len(x) - 1):
        left = (y[i]) / (x[i])
        right = (1 - y[i]) / (1 - x[i])
        angle = calculate_signed_angle(left, right)
        val = max(angle * orientation, 0)
        slopes.append((left, right, val))
        angles.append(val)

    # find the index and value of the max element in ddy_abs
    max_index = np.argmax(angles)
    max_value = angles[max_index]

    # for i, (l, r, v) in enumerate(slopes):
    #     print(f'{i+1}: {l:.2f} {r:.2f} {v:.2f}')
    return max_index + 1, max_value


def find_bends(vals, min_val=0, max_val=1):    
    """
    Finds the knee and elbow in a cumulative distribution of values.

    Parameters
    ----------
    vals : array-like
        Array of values to analyze.
    min_val : float, optional
        Minimum value to consider. Default is 0.
    max_val : float, optional
        Maximum value to consider. Default is 1.

    Returns
    -------
    X : numpy.ndarray
        Array of unique values within the specified range.
    Y : numpy.ndarray
        Cumulative count of samples below each unique value.
    knee_index : int
        Index of the knee.
    knee_value : float
        Value of the knee.
    elbow_index : int
        Index of the elbow.
    elbow_value : float
        Value of the elbow.
    """
    X, Y = count_cum_val(vals, min_val=min_val, max_val=max_val)
    knee_index, knee_value = find_knee(X, Y, knee=True)
    elbow_index, elbow_value = find_knee(X, Y, knee=False)
    return X, Y, knee_index, knee_value, elbow_index, elbow_value


def z_score_filter(vals, min_val=-3, max_val=3):
    """
    Removes values with z-scores below min_val or above max_val.

    Parameters
    ----------
    vals : array-like
        Array of values to filter.
    min_val : float, optional
        Minimum z-score to keep. Default is -3.
    max_val : float, optional
        Maximum z-score to keep. Default is 3.

    Returns
    -------
    numpy.ndarray
        Filtered array of values.
    """
    if min_val is None or max_val is None:
        return vals
    if min_val > max_val:
        raise ValueError('min_val must be less than or equal to max_val')
    zscore = (vals - vals.mean()) / vals.std()
    return vals[(zscore >= min_val) & (zscore <= max_val)]
    