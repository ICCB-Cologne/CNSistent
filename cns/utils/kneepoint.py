import numpy as np
import pandas as pd


# counts how many samples are below each cutoff value
def count_below_lim(vals, min_val=0, max_val=1, steps=1000):
    cutoffs = np.linspace(min_val, max_val, steps + 1)
    #  Finds the indices where elements should be inserted to maintain order, 
    #  effectively counting the number of elements less than or equal to each cutoff.
    counts = np.searchsorted(np.sort(vals), cutoffs, side='right') / len(vals)	
    return cutoffs, counts

# Counts the number of samples below each present values
def count_cum_val(vals, min_val=0, max_val=1):
    vals = np.array(vals)
    vals = vals[(vals >= min_val) & (vals <= max_val)]
    unique_vals, counts = np.unique(vals.astype(np.float32), return_counts=True)
    cumulative_count = np.cumsum(counts) / len(vals)
    return unique_vals, cumulative_count


# concave angles (rising curve) are positive, convex angles (falling curve) are negative
def calculate_signed_angle(s1, s2):
    # Calculate the angle in radians
    angle_radians = np.arctan((s2 - s1) / (1 + s1 * s2))
    
    # Convert the angle to degrees
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees


def get_direction(y):
    if  np.all(np.diff(y) >= 0):
        return 1
    elif np.all(np.diff(y) < 0):
        return -1
    else:
        raise ValueError('Trying to find a knee in a curve that is non-monotonic')


# finds a knee/elbow in the curve when using convex/concave curves
def find_knee(x, y, knee=True):
    if len(x) < 2:
        return -1, np.nan
    y_range = y[-1] - y[0]
    x_range = x[-1] - x[0]
    if y_range == 0:
        return -1, np.nan
    y = (np.array(y) - y[0]) / y_range 
    x = (np.array(x) - x[0]) / x_range
    orientation = get_direction(y)
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


# finds the knee and elbow in a cumulative distribution of values separated by steps between min and max
# vals - values to be counted
# min_val - minimum value to be considered
# max_val - maximum value to be considered
# steps - number of steps between min_val and max_val
def find_bends(vals, min_val=0, max_val=1, steps=1000):    
    X, Y = count_below_lim(vals, min_val=min_val, max_val=max_val, steps=steps)
    knee_index, knee_value = find_knee(X, Y, knee=True)
    elbow_index, elbow_value = find_knee(X, Y, knee=False)
    return X, Y, knee_index, knee_value, elbow_index, elbow_value


# finds the knee and elbow in a cumulative distribution of values separated by steps between min and max
# vals - values to be counted
# min_val - minimum value to be considered
# max_val - maximum value to be considered
# steps - number of steps between min_val and max_val
def find_bends(vals, min_val=0, max_val=1, steps=10000):    
    X, Y = count_below_lim(vals, min_val=min_val, max_val=max_val, steps=steps)
    knee_index, knee_value = find_knee(X, Y, knee=True)
    elbow_index, elbow_value = find_knee(X, Y, knee=False)
    return X, Y, knee_index, knee_value, elbow_index, elbow_value


# removes values with z-scores below min_val or above max_val
def z_score_filter(vals, min_val=-3, max_val=3):
    if min_val is None or max_val is None:
        return vals
    if min_val > max_val:
        raise ValueError('min_val must be less than or equal to max_val')
    zscore = (vals - vals.mean()) / vals.std()
    return vals[(zscore >= min_val) & (zscore <= max_val)]
    