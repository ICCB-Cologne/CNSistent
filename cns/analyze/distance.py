
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
    ids = cns_df["sample_id"].unique()
    chroms = cns_df["chrom"].unique()
    data_df = cns_df.pivot_table(index="sample_id", columns="chrom", values="total_cn", aggfunc=list)
    arr  = np.zeros((len(ids), len(ids)))
    total = len(ids) **2
    for i1, id1 in enumerate(ids):
        for i2, id2 in enumerate(ids[:i1]):
            dists = []
            for chrom in chroms:
                if print_progress:
                    print(f"{(i1*len(ids) + i2)/total:.2%}", end="\r")
                    
                cns1 = data_df.loc[id1, chrom]
                cns2 = data_df.loc[id2, chrom]

                if len(cns1) == 0 or len(cns2) == 0:
                    raise ValueError(f"Wassserstein distance calculation failed on {chrom} for {id1} and {id2}. Both samples must have data to calculate Wasserstein distance.")
                if len(cns1) != len(cns2):
                    raise ValueError(f"Wassserstein distance calculation failed on {chrom} for {id1} and {id2}. Both samples must have the same number of data points to calculate Wasserstein distance.")
                
                dist = calc_wass_distance(cns1, cns2)
                dists.append(dist)
            arr[i1, i2] = np.mean(dists)
            arr[i2, i1] = np.mean(dists)
    return pd.DataFrame(arr, index=ids, columns=ids)
