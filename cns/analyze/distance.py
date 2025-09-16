import pandas as pd
import numpy as np


def manhattan_distance(x, y):
    return np.abs(x - y).sum()


def euclidean_distance(x, y):
    return np.sqrt(((x - y) ** 2).sum())


def wasserstein_distance(x, y):
    x_cumsum = np.cumsum(np.insert(x, 0, 0))
    y_cumsum = np.cumsum(np.insert(y, 0, 0))
    return np.sum(np.abs(x_cumsum - y_cumsum))


def _get_distance_function(dist_type):
    if dist_type == 'manhattan':
        return manhattan_distance
    elif dist_type == 'euclidean':
        return euclidean_distance
    elif dist_type == 'wasserstein':
        return wasserstein_distance
    raise ValueError(f"Unsupported distance type: {dist_type}. Supported types are 'manhattan', 'euclidean', and 'wasserstein'.")


def calc_distances(cns_df, cn_column, dist_type='manhattan'):
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
    dist_type : str, optional
        The type of distance to calculate. One of 'manhattan', 'euclidean', 'wasserstein. Default is 'manhattan'.
    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the pairwise L1 distances between samples, with sample IDs as both index and columns.  
    """
    dist_func = _get_distance_function(dist_type)

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
            dist[i, j] = dist_func(values[i], values[j])
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

