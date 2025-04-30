
import scipy.stats as stats
import pandas as pd
import numpy as np

def calc_wass_distance(sample_A, sample_B, cn_column):
    """
    Calculates the Wasserstein distance between two samples, per chromosomes.

    Parameters
    ----------
    sample_A : pandas.DataFrame
        DataFrame containing CNS data for sample A.
    sample_B : pandas.DataFrame
        DataFrame containing CNS data for sample B.
    cn_column : str
        Column name for copy number data.

    Returns
    -------
    pandas.Series
        Series with the Wasserstein distance for each chromosome.
    """
    
    if sample_A.shape[0] == 0 or sample_B.shape[0] == 0:
        raise ValueError("Both samples must have data to calculate Wasserstein distance.")
    
    # Get common chromosomes
    chroms_A = set(sample_A['chrom'].unique())
    chroms_B = set(sample_B['chrom'].unique())
    common_chroms = chroms_A.intersection(chroms_B)
    
    if not common_chroms:
        raise ValueError("No common chromosomes found between samples")
    
    # Pre-compute groups for better performance
    groups_A = dict(list(sample_A.groupby('chrom')))
    groups_B = dict(list(sample_B.groupby('chrom')))
    
    # Dictionary to store results
    distances = {}
    
    for chrom in common_chroms:
        group_A = groups_A[chrom]
        group_B = groups_B[chrom]
        
        val_count = max(len(group_A), len(group_B))  # Use the larger count for better resolution
        if val_count == 0:
            distances[chrom] = 0
            continue

        vals_A = group_A[cn_column].values
        vals_B = group_B[cn_column].values

        cdf_A = np.concatenate(([0], np.cumsum(vals_A)))
        norm_A = np.sum(vals_A)
        cdf_B = np.concatenate(([0], np.cumsum(vals_B)))
        norm_B = np.sum(vals_B)
        if norm_A == 0 and norm_B == 0:
            cdf_diff = np.zeros_like(cdf_A)
        elif norm_A == 0:
            cdf_diff = cdf_B / norm_B
        elif norm_B == 0:
            cdf_diff = cdf_A / norm_A
        else:
            cdf_diff = np.abs(cdf_A / norm_A - cdf_B / norm_B)
        area_diff = np.trapz(cdf_diff, dx=1/val_count)

        distances[chrom] = area_diff
    
    return pd.Series(distances)