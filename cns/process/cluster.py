import numpy as np
from numba import jit
from cns.utils import hg19

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


# Merge breakpoints that are within a certain distance, extend the first and last to the telomeres
def created_merged_segs(dict_start, dist, assembly=hg19, extend=True):
    res = []
    for chrom, old_breaks in dict_start.items():
        clusters = breaks_to_clusters(old_breaks)
        merged = merge_clusters(clusters, dist)
        new_breaks = clusters_to_breaks(merged)
        if extend and new_breaks[0] < dist:
            new_breaks[0] = 0
        elif new_breaks[0] != 0:
            new_breaks = [0] + new_breaks
        if extend and new_breaks[-1] + dist > assembly.chr_lens[chrom]:
            new_breaks[-1] = assembly.chr_lens[chrom]
        elif new_breaks[-1] != assembly.chr_lens[chrom]:
            new_breaks = new_breaks + [assembly.chr_lens[chrom]]
        
        for i in range(len(new_breaks) - 1):
            res.append((chrom, new_breaks[i] + 1, new_breaks[i + 1]))
    return res