import numpy as np
from cns.process.breakpoints import calc_arm_breaks, calc_cytoband_breaks, create_step_breaks
from cns.utils.conversions import breaks_to_segments, cytobands_to_df, genome_to_segments, tuples_to_segments
from cns.utils.assemblies import hg19
from cns.utils import load_segments


def count_segments(segs):
    return sum([len(chr_segs) for chr_segs in segs.values()])


def do_segments_overlap(segs, is_sorted=False):
    for chr, chr_segs in segs.items():
        if not is_sorted:
            chr_segs.sort(key=lambda x: (x[0]))
        # Check for overlaps
        for i in range(len(chr_segs) - 1):
            current_end = chr_segs[i][1]
            next_start = chr_segs[i + 1][0]
            if current_end > next_start:
                return True
    return False


def find_overlaps(segs, is_sorted=False):
    overlaps = {}
    for chr, chr_segs in segs.items():
        if not is_sorted:
            chr_segs.sort(key=lambda x: (x[0]))
        # Iterate through all pairs of triplets to check for overlap
        n = len(chr_segs)
        for i in range(n):
            current_end = chr_segs[i][1]
            for j in range(i + 1, n):
                next_start = chr_segs[j][0]
                if current_end <= next_start:
                    break
                # Store the overlap along with the group identifiers
                if chr not in overlaps:
                    overlaps[chr] = []
                overlaps[chr].append((next_start, current_end))

    return overlaps


def merge_segments(segs, is_sorted=False):
    merged = {}
    for chr, chr_segs in segs.items():
        if len(chr_segs) == 0:
            merged[chr] = []
            continue

        if not is_sorted:
            chr_segs.sort(key=lambda x: (x[0]))

        merged[chr] = [chr_segs[0]]

        for current in chr_segs[1:]:
            last_start, last_end = merged[chr][-1][0], merged[chr][-1][1]
            last_name = merged[chr][-1][2] if len(merged[chr][-1]) > 2 else None

            # If the current segment starts at the end of the last one
            if current[0] <= last_end:
                merged[chr][-1] = (last_start, current[1]) if last_name is None else (last_start, current[1], last_name)
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
        for seg_a in chr_segs_a:
            new_start = seg_a[0]
            name_a = seg_a[2] if len(seg_a) > 2 else None
            subsections = []
            for seg_b in chr_segs_b:
                # Skip chr_segs_b that are in a different group or before the current segment in chr_segs_a
                if seg_b[1] < new_start:
                    continue
                # Break if the segment in chr_segs_b is beyond the current segment in chr_segs_a
                if seg_b[0] > seg_a[1]:
                    break

                # Calculate the difference if there's an overlap
                if seg_b[0] <= new_start < seg_b[1]:
                    # If chr_segs_a starts within chr_segs_b, move its start to the end of chr_segs_b
                    new_start = seg_b[1]
                elif new_start < seg_b[0] and seg_a[1] > seg_b[0]:
                    # If chr_segs_a overlaps the start of chr_segs_b, add the non-overlapping part to the difference
                    new_seg = (new_start, seg_b[0], name_a) if name_a is not None else (new_start, seg_b[0])
                    subsections.append(new_seg)
                    new_start = seg_b[1]

            # Check if there's any remaining part of chr_segs_a after processing overlaps
            if new_start < seg_a[1]:
                new_seg = (new_start, seg_a[1], name_a) if name_a is not None else (new_start, seg_a[1])
                subsections.append(new_seg)

            if len(subsections) > 1 and len(subsections[0]) > 2:
                subsections = [
                    (subsections[i][0], subsections[i][1], f"{subsections[i][2]}_{i}") for i in range(len(subsections))
                ]

            diffs[chr] += subsections

    return diffs


def filter_min_size(chr_segs, min_size, merge_first=False):
    if merge_first:
        chr_segs = merge_segments(chr_segs)
    return {chr: [seg for seg in chr_segs if seg[1] - seg[0] >= min_size] for chr, chr_segs in chr_segs.items()}


def split_segment(seg_start, seg_end, seg_name, step_size, strategy="scale"):
    length = seg_end - seg_start
    breaks = create_step_breaks(length, step_size, strategy)
    breaks = (np.array(breaks) + seg_start).tolist()
    if seg_name == None:
        res = [(breaks[i], breaks[i + 1]) for i in range(len(breaks) - 1)]
    else:
        res = [(breaks[i], breaks[i + 1], f"{seg_name}_{i}") for i in range(len(breaks) - 1)]
    return res


def split_segments(segments, step_size, strategy="scale"):
    res = {}
    for chr, chr_segs in segments.items():
        res[chr] = []
        for seg in chr_segs:
            seg_name = seg[2] if len(seg) > 2 else None
            res[chr] += split_segment(seg[0], seg[1], seg_name, step_size, strategy)
    return res


def regions_select(select, assembly=hg19):
    if select == "arms":
        arm_breaks = calc_arm_breaks(assembly)
        return {
            chrom: [(breaks[0], breaks[1], f"{chrom}p"), (breaks[1], breaks[2], f"{chrom}q")]
            for chrom, breaks in arm_breaks.items()
        }
    elif select == "bands":
        bands_df = cytobands_to_df(assembly.cytobands)
        return {
            chrom: [(start, end, name) for start, end, name in zip(subset["start"], subset["end"], subset["name"])]
            for chrom, subset in bands_df.groupby("chrom")
        }
    elif select == "":
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
        res = filter_min_size(res, filter_size, False)
    if remove != None:
        if filter_size > 0:
            remove = filter_min_size(remove, filter_size, True)
        res = segment_difference(res, remove)
        if filter_size > 0:
            res = filter_min_size(res, filter_size, False)
    return res