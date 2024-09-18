import numpy as np
from cns.process.breakpoints import calc_arm_breaks, calc_cytoband_breaks, create_step_breaks
from cns.utils.assemblies import hg19
from cns.utils.conversions import breaks_to_segments, genome_to_segments, tuples_to_segments
from cns.utils.files import load_segments


def do_segments_overlap(segs, sorted=False):
    for chr, chr_segs in segs.items():
        if not sorted:
            chr_segs.sort(key=lambda x: (x[0]))    
        # Check for overlaps
        for i in range(len(chr_segs) - 1):
            current_end = chr_segs[i][1]
            next_start = chr_segs[i+1][0]
            if  current_end > next_start:
                return True
    return False


def find_overlaps(segs, sorted=False):
    overlaps = {} 
    for chr, chr_segs in segs.items():
        if not sorted:
            chr_segs.sort(key=lambda x: (x[0]))    
        # Iterate through all pairs of triplets to check for overlap
        n = len(chr_segs)
        for i in range(n):
            current_end = chr_segs[i][1]
            for j in range(i+1, n):
                next_start = chr_segs[j][0]
                if current_end <= next_start:
                    break                
                # Store the overlap along with the group identifiers
                if chr not in overlaps:
                    overlaps[chr] = []
                overlaps[chr].append((next_start, current_end))
    
    return overlaps


def merge_segments(segs, sort=False):
    merged = {}
    for chr, chr_segs in segs.items():
        if not sorted:
            chr_segs.sort(key=lambda x: (x[0]))   

        merged[chr] = [chr_segs[0]]

        for current in chr_segs[1:]:
            last_start, last_end = merged[chr][-1][0], merged[chr][-1][1]

            # If the current segment starts at the end of the last one plus one
            if current[0] <= last_end:
                merged[chr][-1] = (last_start, current[1])
            else:
                # Add the current segment as is
                merged[chr].append(current)

    return merged


def segment_union(segs_a, segs_b):
    # Combine and sort the segments first by group, then by start time
    keys = set(segs_a.keys()).union(set(segs_b.keys()))
    new_segs = {}
    for key in keys:
        new_segs[key] = segs_a.get(key, []) + segs_b.get(key, [])
        new_segs[key].sort(key=lambda x: (x[0]))
    merged = merge_segments(new_segs)
    return merged


def segment_difference(segs_a, segs_b, sorted=False):
    diffs = {}  
    for chr, chr_segs_a in segs_a.items():
        if chr not in segs_b:
            diffs[chr] = chr_segs_a
            continue
        chr_segs_b = segs_b[chr]

        if not sorted:
            chr_segs_a.sort(key=lambda x: (x[0]))
            chr_segs_b.sort(key=lambda x: (x[0]))    

        diffs[chr] = []

        # Iterate through each segment in chr_segs_a
        for start_a, end_a in chr_segs_a:
            new_start = start_a
            for start_b, end_b in chr_segs_b:                
                # Skip chr_segs_b that are in a different group or before the current segment in chr_segs_a
                if end_b < new_start:
                    continue
                # Break if the segment in chr_segs_b is beyond the current segment in chr_segs_a
                if  start_b > end_a:
                    break
                
                # Calculate the difference if there's an overlap
                if start_b <= new_start < end_b:
                    # If chr_segs_a starts within chr_segs_b, move its start to the end of chr_segs_b
                    new_start = end_b
                elif new_start < start_b and end_a > start_b:
                    # If chr_segs_a overlaps the start of chr_segs_b, add the non-overlapping part to the difference
                    diffs[chr].append((new_start, start_b))
                    new_start = end_b
            
            
            # Check if there's any remaining part of chr_segs_a after processing overlaps
            if new_start < end_a:
                diffs[chr].append((new_start, end_a))
            # Reset index_b for the next iteration through chr_segs_a

    return diffs


def filter_min_size(chr_segs, min_size):
    return { chr : [seg for seg in chr_segs if seg[1] - seg[0] >= min_size] for chr, chr_segs in chr_segs.items() }


def split_segment(seg_start, seg_end, step_size, strategy="scale"):
    length = seg_end - seg_start
    breaks = create_step_breaks(length, step_size, strategy)
    breaks = (np.array(breaks) + seg_start).tolist()
    res = [(breaks[i], breaks[i + 1]) for i in range(len(breaks) - 1)]
    return res


def split_segments(segments, step_size, strategy="scale"):
    res = {}
    for chr, chr_segs in segments.items():
        res[chr] = []
        for seg_start, seg_end in chr_segs:
            res[chr] += split_segment(seg_start, seg_end , step_size, strategy)
    return res


def regions_select(select, assembly=hg19):
    if select == "arms":
        breaks = calc_arm_breaks(assembly)
        return breaks_to_segments(breaks)
    elif select == "bands":
        breaks = calc_cytoband_breaks(assembly)
        return breaks_to_segments(breaks)
    elif select =="":
        return genome_to_segments(assembly)
    else:
        return load_segments(select)


def regions_remove(remove, assembly=hg19):
    if remove == "gaps":
        breaks = assembly.gaps
        return tuples_to_segments(breaks)
    elif remove == "":
        return None
    else:
        return load_segments(remove)


def get_genome_segments(select, remove=None, filter_size=0):
    res = select
    if filter_size > 0:
        res = filter_min_size(res, filter_size)
    if remove != None:
        if filter_size > 0:
            remove = filter_min_size(remove, filter_size)
        res = segment_difference(res, remove)
        if filter_size > 0:
            res = filter_min_size(res, filter_size)
    return res



def calc_chr_sizes(chr_segs_df):
    return { chr : sum([seg[1] - seg[0] for seg in chr_segs]) for chr, chr_segs in chr_segs_df.items() }


def calc_group_sizes(segs_df, assembly=hg19):
    chr_sizes = calc_chr_sizes(segs_df)
    groups = { "sex" : assembly.sex_names, "aut" : assembly.aut_names, "all" : assembly.chr_names }
    return { key : sum([chr_sizes[chrom] if chrom in chr_sizes else 0 for chrom in group]) for key, group in groups.items() } 