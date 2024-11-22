import numpy as np
from numba import jit
from cns.utils.conversions import segments_to_breaks
from cns.utils.logging import log_info

@jit(nopython=True)
def _merge_clusters(clusters, threshold):
    merged = True
    while merged:
        merged = False
        i = 0
        while i < len(clusters)-1:
            # Calculate the centroid of the current cluster
            centroid_i = clusters[i, 0]
            if centroid_i < 0:
                i += 1
                continue

            # Look ahead to merge with any close clusters
            for j in range(i+1, len(clusters)):
                # Calculate the centroid of the next cluster
                centroid_j = clusters[j, 0]
                if centroid_j < 0:
                    continue
                
                # Check if the distance between centroids is below the threshold
                if abs(centroid_i - centroid_j) < threshold:
                    # Merge clusters i and j
                    new_count = clusters[i, 1] + clusters[j, 1]
                    new_center = (clusters[i, 0] * clusters[i, 1] + clusters[j, 0] * clusters[j, 1]) / new_count
                    clusters[i, 0] = new_center
                    clusters[i, 1] = new_count
                    # Remove the merged cluster
                    clusters[j, 0] = -1
                    clusters[j, 1] = -1
                    merged = True
                    # After merging, break to restart the process (to consider new distances)
                    break
                
            i += 1

    filtered = clusters[clusters[:, 0] >= 0]
    rounded = np.round(filtered).astype(np.int64)
    return rounded


def _breaks_to_clusters(chrom_breaks):
    return np.array([[val, 1] for val in chrom_breaks], dtype=np.float64)


def _break_count(breaks):
    return sum(len(values) for values in breaks.values())


def _insert_breaks_into_segments(segs, breaks, min_dist=0):
    res = {}
    for chrom, chrom_segs in segs.items():
        if chrom not in breaks:
            res[chrom] = chrom_segs
            continue
        res[chrom] = []
        for seg in chrom_segs:
            seg_breaks = []
            start = seg[0]
            end = seg[1]
            name = seg[2] if len(seg) > 2 else None
            for br in breaks[chrom]:
                if start + min_dist <= br < end - min_dist:
                    seg_breaks.append(br)
            if len(seg_breaks) == 0:
                new_seg = (start, end, name) if name != None else (start, end)
                res[chrom].append(new_seg)
            else:
                seg_breaks = [start] + sorted(set(seg_breaks)) + [end]
                for i in range(len(seg_breaks) - 1):
                    if seg_breaks[i + 1] - seg_breaks[i] > 0:
                        new_seg = (seg_breaks[i], seg_breaks[i + 1], name) if name != None else (seg_breaks[i], seg_breaks[i + 1])
                        res[chrom].append(new_seg)
    return res


def _get_breaks_inside_segments(segs, breaks, min_dist = 0):
    res = { chrom: [] for chrom in breaks }
    for chrom, chrom_segs in segs.items():
        if chrom not in breaks:
            continue
        for seg in chrom_segs:
            for br in breaks[chrom]:
                if seg[0] + min_dist <= br < seg[1] - min_dist:
                    res[chrom].append(br)
    return res


def clusters_to_breaks(clusters):
    chrom_breaks = []
    for value in clusters:
        chrom_breaks.append(value[0])
    return chrom_breaks


def calc_clusters(dict_start, dist):
    res = {}
    for chrom, old_breaks in dict_start.items():
        new_breaks = []
        if len(old_breaks) > 0:
            clusters = _breaks_to_clusters(old_breaks)
            merged = _merge_clusters(clusters, dist) if dist > 0 else clusters
            new_breaks = clusters_to_breaks(merged)
        res[chrom] = new_breaks
    return res


def cluster_segments(segments, clust_dist, print_info=False):
    breaks = segments_to_breaks(segments)
    return cluster_breaks(breaks, clust_dist, print_info)


def cluster_breaks(breaks, clust_dist, keep_ends = True, print_info=False):    
    break_count = _break_count(breaks)
    log_info(print_info, f"Reducing {break_count} breakpoints with merge distance {clust_dist} ... ")

    res = calc_clusters(breaks, clust_dist)
    if keep_ends:
        for chrom, chrom_breaks in breaks.items():
            if res[chrom][0] != chrom_breaks[0]:
                res[chrom].insert(0, chrom_breaks[0])
            if res[chrom][-1] != chrom_breaks[-1]:
                res[chrom].append(chrom_breaks[-1])
    new_count = _break_count(res)
    log_info(print_info, f"Removed to {break_count - new_count} breakpoints by distance merge.")
    log_info(print_info, f"Resulting breakpoints: {new_count}")

    return res
