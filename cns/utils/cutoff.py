import numpy as np
import pandas as pd


# counts how many samples are below each cutoff value
def count_below_lim(vals, min_val=0, max_val=1, steps=1000):
    cutoffs = np.linspace(min_val, max_val, steps + 1)
    #  Finds the indices where elements should be inserted to maintain order, 
    #  effectively counting the number of elements less than or equal to each cutoff.
    counts = np.searchsorted(np.sort(vals), cutoffs, side='right') / len(vals)	
    return cutoffs, counts


# concave angles (rising curve) are positive, convex angles (falling curve) are negative
def calculate_signed_angle(s1, s2):
    # Calculate the angle in radians
    angle_radians = np.arctan((s2 - s1) / (1 + s1 * s2))
    
    # Convert the angle to degrees
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees


# finds a knee/elbow in the curve when using convex/concave curves
def find_knee(x, y, convex=True):
    y_range = y[-1] - y[0]
    x_range = x[-1] - x[0]
    norm_factor = x_range / y_range

    # calculate the difference between slopes on the left and right side of each point 
    angles = []
    slopes = []
    for i in range(1, len(x) - 1):
        left = (y[i] - y[0]) / (x[i] - x[0]) * norm_factor
        right = (y[-1] - y[i]) / (x[-1] - x[i]) * norm_factor
        angle = calculate_signed_angle(left, right)
        val = max(angle * -1, 0) if convex else max(angle, 0)
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
# steps - number of cutoffs between min and max, such that the distribution is counts the vals <= cutoff
# dist - number of points left and right of a cutoff to be considered when calculating the slope
# allow_pad - if True, the slope is calculated for all points, 
#   otherwise the slope is calculated only for points that are not within dist from the beginning or end
def find_bends(vals, min_val=0, max_val=1, steps=1000, dist=10, allow_boundary=True):    
    cutoffs, counts, delta_x = count_below_lim(vals, min_val=min_val, max_val=max_val, steps=steps)
    if delta_x <= 0:
        return cutoffs, counts, -1, min_val, -1, min_val
    knee_index, knee_value = find_knee(counts, delta_x, convex=True, dist=dist, allow_boundary=allow_boundary)
    elbow_index, elbow_value = find_knee(counts, delta_x, convex=False, dist=dist, allow_boundary=allow_boundary)
    return cutoffs, counts, knee_index, knee_value, elbow_index, elbow_value