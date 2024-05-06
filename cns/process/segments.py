

import numpy as np
from cns.process.breakpoints import split_segments
from cns.utils.assemblies import hg19


def genome_to_segments(chr_lens):
    regions = []
    for chrom, len in chr_lens.items():
        regions.append((chrom, 0, len))
    return regions


def breaks_to_segments(breakpoints):
    segments = []
    for chrom, breaks in breakpoints:
        last_break = len(breaks) - 1
        for i in range(last_break):
            segments.append((chrom, breaks[i], breaks[i + 1]))
    return segments


def regions_to_segments(regions, change_coords = False):
    segments = []
    for chrom, start, end in regions[["chrom", "start", "end"]].values:
        segments.append((chrom, start - 1 if change_coords else start, end))
    return segments


def tuples_to_segments(tuples):
    segments = []
    if len(tuples) > 0 and len(tuples[0]) >= 3:
        for tuple in tuples:
            segments.append((tuple[0], tuple[1], tuple[2]))
    return segments


def do_segments_overlap(segs, sorted=False):
    # Sort segments by group, then by start time
    if not sorted:
        segs.sort(key=lambda x: (x[0], x[1]))    
    # Check for overlaps
    for i in range(len(segs) - 1):
        current_group, current_start, current_end = segs[i]
        next_group, next_start, next_end = segs[i+1]
        
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
        group1, start1, end1 = segs[i]
        for j in range(i+1, n):
            group2, start2, end2 = segs[j]
            if group1 != group2 or end1 <= start2:
                break
            
            # Store the overlap along with the group identifiers
            overlaps.append((group1, start2, end1))
    
    return overlaps


def merge_segments(segs):
    # Sort segments by start time
    segs.sort(key=lambda x: (x[0], x[1]))

    merged = [segs[0]]

    for current in segs[1:]:
        last_group, last_start, last_end = merged[-1]

        # If the current segment starts at the end of the last one plus one
        if current[1] <= last_end and current[0] == last_group:
            # Merge the two segments
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


def get_genome_segments(select, bin_size=0, remove=None, filter_size=0):
    res = select
    if filter_size > 0:
        res = filter_min_size(res, filter_size)
    if remove != None:
        if filter_size > 0:
            remove = filter_min_size(remove, filter_size)
        res = segment_difference(res, remove)
        if filter_size > 0:
            res = filter_min_size(res, filter_size)
    if bin_size > 0:
        res = split_segments(res, bin_size)
    return res


def add_seg_info(cns_df, assembly=hg19):
    cns_df = cns_df.copy()
    cns_df["length"] = (cns_df["end"] - cns_df["start"]).astype(np.uint32)
    cns_df["mid"] = cns_df["start"] + cns_df["length"] // 2
    cns_df["cum_mid"] = cns_df["mid"] + cns_df.apply(lambda x: assembly.chr_starts[x["chrom"]], axis=1)
    if "major_cn" and "minor_cn" in cns_df:
        cns_df["total_cn"] = cns_df["major_cn"] + cns_df["minor_cn"]
    if "cn_a" and "cn_b" in cns_df:
        cns_df["total_cn"] = cns_df["cn_a"] + cns_df["cn_b"]
    # order by cum_mid
    cns_df = cns_df.sort_values(by=["sample_id", "cum_mid"])
    return cns_df
