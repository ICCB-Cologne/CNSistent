import numpy as np
from numba import jit
from cns.process.segments import get_consecutive_segs
from cns.utils.conversions import values_count
from cns.utils.logging import log_info

@jit(nopython=True)
def _merge_clusters(clusters, threshold):
    i = 0
    while i < len(clusters)-1:
        # Calculate the centroid of the current cluster
        start_i = clusters[i, 2]
        if start_i < 0:
            i += 1
            continue

        # Look ahead to merge with any close clusters
        for j in range(i+1, len(clusters)):
            # Calculate the centroid of the next cluster
            start_j = clusters[j, 2]
            if start_j < 0:
                continue
            
            # Check if the distance between centroids is below the threshold
            if abs(start_i - start_j) <= threshold:
                # Merge clusters i and j
                new_count = clusters[i, 1] + clusters[j, 1]
                new_center = (clusters[i, 0] * clusters[i, 1] + clusters[j, 0] * clusters[j, 1]) / new_count
                clusters[i, 0] = new_center
                clusters[i, 1] = new_count
                # Remove the merged cluster
                clusters[j, :] = -1
            
        i += 1

    filtered = clusters[clusters[:, 0] >= 0]
    rounded = np.round(filtered).astype(np.int64)
    return rounded


def _breaks_to_clusters(chrom_breaks):
    return np.array([[val, 1, val] for val in chrom_breaks], dtype=np.float64)


def _clusters_to_breaks(clusters):
    chrom_breaks = []
    for value in clusters:
        chrom_breaks.append(value[0])
    return chrom_breaks


def _cluster_breaks_list(breaks, clust_dist, keep_ends = True, print_info=False):
    if len(breaks) < 2:
        return breaks
    if len(breaks) == 2 and keep_ends:
        return breaks
    
    if keep_ends:
        clusters = _breaks_to_clusters(breaks[1:-1])	
    else:
        clusters = _breaks_to_clusters(breaks)
    merged = _merge_clusters(clusters, clust_dist)
    new_breaks = _clusters_to_breaks(merged)
    if keep_ends:
        new_breaks = [breaks[0]] + new_breaks + [breaks[-1]]
    return new_breaks


def _get_break_list(segs):
    if len(segs) == 0:
        return []
    res = [segs[0][0]]
    for seg in segs:
        res.append(seg[1])
    return res


def _extend_segs(chrom, breaks, chr_segs):
    offset = len(chr_segs)
    for i in range(len(breaks) - 1):
        chr_segs.append((breaks[i], breaks[i + 1], f'{chrom}_{offset+i}'))


def cluster_segments(input_segs, clust_dist, keep_ends = True, print_info=False):
    """
    Clusters segments based on the specified clustering distance.

    Parameters
    ----------
    input_segs : dict
        Dictionary of input segments with chromosome names as keys and list of segments as values.
    clust_dist : int
        Clustering distance to use for merging segments.
    keep_ends : bool, optional
        If True, keeps the ends of the segments. Default is True.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    dict
        Dictionary of clustered segments with chromosome names as keys and list of segments as values.
    """
    seg_count = values_count(input_segs)
    log_info(print_info, f"Merging {seg_count} segments with merge distance {clust_dist} ... ")
    res = {}
    for chrom, chrom_segs in input_segs.items():
        cons = get_consecutive_segs(chrom_segs)        
        chr_segs = []
        for neigbours in cons:
            if clust_dist > neigbours[-1][1] - neigbours[0][0]:
                chr_segs.append((neigbours[0][0], neigbours[-1][1], f'{chrom}_{len(chr_segs)}'))
            else:
                seg_breaks = _get_break_list(neigbours)
                breaks = _cluster_breaks_list(seg_breaks, clust_dist, keep_ends, print_info)
                _extend_segs(chrom, breaks, chr_segs)       
        res[chrom] = chr_segs

    new_count = values_count(res)
    log_info(print_info, f"Removed {seg_count - new_count} segments by distance merge.")
    log_info(print_info, f"Resulting segment count: {new_count}")
    return res