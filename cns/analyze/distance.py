
import pandas as pd
import numpy as np

def calc_wass_distance(vals_A, vals_B, norm_vals=True, norm_count=True):
    """
    Calculates the Wasserstein distance between two sequences of the same lenght

    Parameters
    ----------
    vals_A : array-like
        First sequence of values.
    vals_B : array-like
        Second sequence of values.

    Returns
    -------
    float
        The Wasserstein distance between the two sequences.
    """
    count_norm = 1/len(vals_A) if norm_count else 1

    cdf_A = np.concatenate(([0], np.cumsum(vals_A)))
    norm_A = np.sum(vals_A) if norm_vals else 1
    cdf_B = np.concatenate(([0], np.cumsum(vals_B)))
    norm_B = np.sum(vals_B) if norm_vals else 1
    if norm_A == 0 and norm_B == 0:
        cdf_diff = np.zeros_like(cdf_A)
    elif norm_A == 0:
        cdf_diff = cdf_B / norm_B
    elif norm_B == 0:
        cdf_diff = cdf_A / norm_A
    else:
        cdf_diff = np.abs(cdf_A / norm_A - cdf_B / norm_B)
    area_diff = np.trapz(cdf_diff, dx=count_norm)

    return area_diff

    
def calc_distances(cns_df, print_progress=False):
    ids = np.array(sorted(cns_df["sample_id"].unique()))
    chroms = np.array(sorted(cns_df["chrom"].unique()))
    n_ids = len(ids)
    arr = np.zeros((n_ids, n_ids))
    total = n_ids ** 2

    # Build a 2D numpy array for each chromosome: rows=samples, columns=CN bins
    chrom_arrays = {}
    for chrom in chroms:
        chrom_df = cns_df[cns_df["chrom"] == chrom]
        # Optionally sort by 'start' if available for consistent bin order
        chrom_df = chrom_df.sort_values(by=['sample_id', 'start'])
        # Pivot to (sample_id, bin) -> total_cn
        pivot = chrom_df.pivot(index="sample_id", columns="start", values="total_cn")
        # Ensure all samples are present and order matches ids
        pivot = pivot.reindex(ids)
        # Each row is a CN vector for a sample for this chromosome
        chrom_arrays[chrom] = np.vstack(pivot.values)

    for i1, id1 in enumerate(ids):
        for i2, id2 in enumerate(ids[:i1]):
            dists = []
            for chrom in chroms:
                cns1 = chrom_arrays[chrom][i1]
                cns2 = chrom_arrays[chrom][i2]
                if len(cns1) == 0 or len(cns2) == 0:
                    raise ValueError(f"Wasserstein distance calculation failed on {chrom} for {id1} and {id2}. Both samples must have data.")
                if len(cns1) != len(cns2):
                    raise ValueError(f"Wasserstein distance calculation failed on {chrom} for {id1} and {id2}. Both samples must have the same number of data points.")
                dist = calc_wass_distance(cns1, cns2)
                dists.append(dist)
            arr[i1, i2] = np.mean(dists)
            arr[i2, i1] = arr[i1, i2]
            if print_progress:
                print(f"{(i1*n_ids + i2)/total:.2%}", end="\r")
    return pd.DataFrame(arr, index=ids, columns=ids)
