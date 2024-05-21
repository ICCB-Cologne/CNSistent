import numpy as np
from cns.process.breakpoints import create_step_breaks


def do_segments_overlap(segs, sorted=False):
    # Sort segments by group, then by start time
    if not sorted:
        segs.sort(key=lambda x: (x[0], x[1]))    
    # Check for overlaps
    for i in range(len(segs) - 1):
        current_group, current_end = segs[i][0], segs[i][2]
        next_group, next_start = segs[i+1][0], segs[i+1][1]
        
        # Check if they are in the same group and overlap
        if current_group == next_group and current_end > next_start:
            return True
    return False


def find_overlaps(segs, sorted=False):
    if not sorted:
        segs.sort(key=lambda x: (x[0], x[1]))    
    overlaps = []    
    # Iterate through all pairs of triplets to check for overlap
    n = len(segs)
    for i in range(n):
        current_group, current_end = segs[i][0], segs[i][2]
        for j in range(i+1, n):
            next_group, next_start = segs[j][0], segs[j][1]
            if current_group != next_group or current_end <= next_start:
                break
            
            # Store the overlap along with the group identifiers
            overlaps.append((current_group, next_start, current_end))
    
    return overlaps


def merge_segments(segs):
    # Sort segments by start time
    segs.sort(key=lambda x: (x[0], x[1]))

    merged = [segs[0]]

    for current in segs[1:]:
        last_group, last_start, last_end = merged[-1][0], merged[-1][1], merged[-1][2]

        # If the current segment starts at the end of the last one plus one
        if current[1] <= last_end and current[0] == last_group:
            merged[-1] = (last_group, last_start, current[2])
        else:
            # Add the current segment as is
            merged.append(current)

    return merged


def segment_union(segs_a, segs_b):
    # Combine and sort the segments first by group, then by start time
    segs = segs_a + segs_b
    merged = merge_segments(segs)
    return merged


def segment_difference(segs_a, segs_b, sorted=False):
    differences = []

    if not sorted:
        segs_a.sort(key=lambda x: (x[0], x[1]))
        segs_b.sort(key=lambda x: (x[0], x[1]))    
    
    index_b = 0
    # Iterate through each segment in segs_a
    for group_a, start_a, end_a in segs_a:
        new_start = start_a
        while index_b < len(segs_b):
            group_b, start_b, end_b = segs_b[index_b]
            
            # Skip segs_b that are in a different group or before the current segment in segs_a
            if group_b < group_a or end_b < new_start:
                index_b += 1
                continue
            # Break if the segment in segs_b is beyond the current segment in segs_a
            if group_b > group_a or start_b > end_a:
                break
            
            # Calculate the difference if there's an overlap
            if start_b <= new_start < end_b:
                # If segs_a starts within segs_b, move its start to the end of segs_b
                new_start = end_b
            elif new_start < start_b and end_a > start_b:
                # If segs_a overlaps the start of segs_b, add the non-overlapping part to the difference
                differences.append((group_a, new_start, start_b))
                new_start = end_b
            
            index_b += 1
        
        # Check if there's any remaining part of segs_a after processing overlaps
        if new_start < end_a:
            differences.append((group_a, new_start, end_a))
        # Reset index_b for the next iteration through segs_a
        index_b = 0

    return differences


def filter_min_size(segs, min_size):
    return [seg for seg in segs if seg[2] - seg[1] >= min_size]


def split_segment(segment, step_size, equidisant=True):
    chrom_no, seg_start, seg_end = segment
    length = seg_end - seg_start
    breaks = create_step_breaks(length, step_size, equidisant)
    breaks = (np.array(breaks) + seg_start).tolist()
    res = [(chrom_no, breaks[i], breaks[i + 1]) for i in range(len(breaks) - 1)]
    return res


def split_segments(segments, step_size, equidisant=True):
    res = []
    for segment in segments:
        res += split_segment(segment, step_size, equidisant)
    return res