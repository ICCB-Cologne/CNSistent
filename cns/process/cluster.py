import numpy as np
from numba import jit
from cns.process.breakpoints import get_breaks, get_breaks_in_segments, insert_breaks_in_segments
from cns.utils import hg19
from cns.utils.logging import log_info

@jit(nopython=True)
def merge_clusters(clusters, threshold):
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


def breaks_to_clusters(chrom_breaks):
    return np.array([[val, 1] for val in chrom_breaks], dtype=np.float64)


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
            clusters = breaks_to_clusters(old_breaks)
            merged = merge_clusters(clusters, dist)
            new_breaks = clusters_to_breaks(merged)
        res[chrom] = new_breaks
    return res


def _break_count(breaks):
    return sum(len(values) for values in breaks.values())


def cluster_within_segments(existing_breaks, segments, clust_dist, print_info=False):    
    break_count = _break_count(existing_breaks)
    log_info(print_info, f"Reducing {break_count} breakpoints with merge distance {clust_dist} ... ")

    existing_breaks = get_breaks_in_segments(segments, existing_breaks)
    new_count = _break_count(existing_breaks)
    if new_count != break_count:
        log_info(print_info, f"Removed {break_count - new_count} outside of segments.")
        break_count = new_count

    res = calc_clusters(existing_breaks, clust_dist)
    new_count = _break_count(existing_breaks)
    log_info(print_info, f"Removed to {break_count - new_count} breakpoints by distance merge.")

    res = insert_breaks_in_segments(segments, res, clust_dist // 2 - 1)

    return res
