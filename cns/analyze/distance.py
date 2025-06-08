import pandas as pd
import numpy as np


def _calc_wass_distances(cns_df, cn_column, norm_length=True, print_progress=False):
    if cn_column not in cns_df.columns:
        raise ValueError(f"Column '{cn_column}' not found in DataFrame.")

    ids = np.array(sorted(cns_df["sample_id"].unique()))
    n_ids = len(ids)
    chroms = np.array(sorted(cns_df["chrom"].unique()))
    n_chroms = len(chroms)

    # Build a 2D numpy array for each chromosome: rows=samples, columns=cumsum CN bins (already normalized)
    if print_progress:
        print("Building 2D arrays for each chromosome...", end="\r")
    chrom_arrays = {}
    cns_df = cns_df.sort_values(by=["sample_id", "chrom", "start"])
    for chrom in chroms:
        chrom_df = cns_df[cns_df["chrom"] == chrom]
        pivot = chrom_df.pivot(index="sample_id", columns="start", values=cn_column)
        cumsum_arr = np.cumsum(np.nan_to_num(pivot.values, nan=0), axis=1)  # each segment "adds" itself
        cumsum_arr = np.concatenate([np.zeros((cumsum_arr.shape[0], 1)), cumsum_arr], axis=1)  # lead with zeros
        # normalize
        row_sums = cumsum_arr[:, -1]
        row_sums[row_sums == 0] = 1  # avoid division by zero, will retain zero in location
        chrom_arrays[chrom] = cumsum_arr / row_sums[:, None]

    # For each chromosome, compute its n_ids x n_ids distance matrix
    if print_progress:
        print("Calculating distance matrices..." + " " * 20)
    chrom_dist_matrices = {}
    for chri, chrom in enumerate(chroms):
        mat = np.zeros((n_ids, n_ids))
        arrs = chrom_arrays[chrom]
        count_norm = 1 / (chrom_arrays[chrom].shape[1]) if norm_length else 1
        step = chri / n_chroms
        for i1 in range(n_ids):
            if print_progress:
                print(f"{i1/(n_ids*n_chroms) + step:.2%}", end="\r")
            for i2 in range(i1 + 1, n_ids):
                cdf_diff = np.abs(arrs[i1] - arrs[i2])
                area_diff = np.trapz(cdf_diff, dx=count_norm)
                mat[i1, i2] = area_diff
                mat[i2, i1] = area_diff  # symmetric matrix
        chrom_dist_matrices[chrom] = mat

    # Calculate the mean distance matrix across chromosomes
    arr = np.mean(np.array(list(chrom_dist_matrices.values())), axis=0)

    return pd.DataFrame(arr, index=ids, columns=ids)


def calc_distances(cns_df, cn_column):
    """
    Calculate the pairwise L1 (Manhattan) distance matrix between samples based on a specified column.
    Each sample is represented as a vector of values for the specified column, with regions as features.
    The values are normalized to proportions by dividing each sample's vector by its sum.
    The function computes the pairwise L1 distances between all samples and returns the result as a DataFrame.
    Parameters
    ----------
    cns_df : pandas.DataFrame
        Input DataFrame containing at least 'sample_id', 'name', and the specified column.
    cn_column : str
        The name of the column in `cns_df` to use for distance calculation.
    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the pairwise L1 distances between samples, with sample IDs as both index and columns.  
    """
    pivot = cns_df.pivot(index="sample_id", columns="name", values=cn_column).fillna(0)
    sample_ids = pivot.index
    values = pivot.values
    # divide by the sum of each row to get proportions
    values = values / values.sum(axis=1, keepdims=True)

    # Calculate pairwise L1 (Manhattan) distances
    n = len(sample_ids)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dist[i, j] = np.abs(values[i] - values[j]).sum()
    # Return as a DataFrame for easier downstream use
    return pd.DataFrame(dist, index=sample_ids, columns=sample_ids)


def calc_chrom_distances(cns_df, cn_column):
    """
    Calculate the pairwise L1 (Manhattan) distance matrix between two samples based on a specified column.
    Provide a value for each chromosome.
       Parameters
    ----------
    cns_df : pandas.DataFrame
        Input DataFrame containing at least 'sample_id', 'name', and the specified column.
    cn_column : str
        The name of the column in `cns_df` to use for distance calculation.
    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the pairwise L1 distances between samples, with sample IDs as both index and columns.  
    """
    if cns_df.sample_id.nunique() != 2:
        raise ValueError("This function only works for two samples.")
    groups = cns_df.groupby("chrom")
    res = {}
    for chrom, group in groups:
        res[chrom] = calc_distances(group, cn_column).iloc[0,1]
    return pd.Series(res)