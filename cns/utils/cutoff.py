import numpy as np
import pandas as pd


# counts how many samples are below each cutoff value
def count_below_lim(vals, min_val=0, max_val=1, steps=1000):
    cutoffs = np.linspace(min_val, max_val, steps + 1)
    #  Finds the indices where elements should be inserted to maintain order, 
    #  effectively counting the number of elements less than or equal to each cutoff.
    counts = np.searchsorted(np.sort(vals), cutoffs, side='right') / len(vals)	
    max_count = np.max(counts)
    delta_x=(max_val - min_val)/(len(cutoffs) * np.max(counts)) if max_count > 0 else 0
    return cutoffs, counts, delta_x # x, y


# concave angles (rising curve) are positive, convex angles (falling curve) are negative
def calculate_signed_angle(s1, s2):
    # Calculate the angle in radians
    angle_radians = np.arctan((s2 - s1) / (1 + s1 * s2))
    
    # Convert the angle to degrees
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees


# finds a knee/elbow in the curve when using covex/concave cureves
def find_knee(vals, delta_x = None, convex=True, dist=1, allow_pad=True):
    if delta_x is None:
        delta_x = 1 / len(vals)

    steps = len(vals) - 1
    y = vals
    for i in range(dist):
        # insert fist value at the beginning of the array dist times
        y = np.insert(y, 0, y[0])
        # insert last value at the end of the array dist times
        y = np.append(y, y[-1])
    
    dy = np.diff(y) / delta_x

    # calculate the difference between slopes on the left and right side of each point 
    ddy_abs = []
    for i in range(steps):
        if not allow_pad and (i < dist or i > steps - dist):
            ddy_abs.append(0)
            continue
        left = np.mean(dy[i:i+dist])
        right = np.mean(dy[i+dist:i+dist*2])
        angle = calculate_signed_angle(right, left)
        val = max(angle, 0) if convex else max(angle * -1, 0)
        ddy_abs.append(val)

    # find the index and value of the max element in ddy_abs
    max_index = np.argmax(ddy_abs)
    max_value = ddy_abs[max_index]

    # print("y", y)
    # print("dy", dy)
    # print("ddy_abs", ddy_abs)

    return max_index, max_value


# finds the knee and elbow in a cumulative distribution of values separated by steps between min and max
# vals - values to be counted
# min_val - minimum value to be considered
# max_val - maximum value to be considered
# steps - number of cutoffs between min and max, such that the distribution is counts the vals <= cutoff
# dist - number of points left and right of a cutoff to be considered when calculating the slope
# allow_pad - if True, the slope is calculated for all points, 
#   otherwise the slope is calculated only for points that are not within dist from the beginning or end
def find_bends(vals, min_val=0, max_val=1, steps=1000, dist=10, allow_pad=True):    
    cutoffs, counts, delta_x = count_below_lim(vals, min_val=min_val, max_val=max_val, steps=steps)
    knee_index, knee_value = find_knee(counts, delta_x, convex=True, dist=dist, allow_pad=allow_pad)
    elbow_index, elbow_value = find_knee(counts, delta_x, convex=False, dist=dist, allow_pad=allow_pad)
    return cutoffs, counts, knee_index, elbow_index