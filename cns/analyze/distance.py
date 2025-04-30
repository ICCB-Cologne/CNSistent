
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
    
    if sample_A.shape[1] != sample_B.shape[1]:
        raise ValueError("Both samples must have the same number of columns.")
    
    groups_A = sample_A.groupby('chrom')
    groups_B = sample_B.groupby('chrom')
    
    # Dictionary to store results
    distances = 0
    
    for chrom, group_A in groups_A:
        # Filter data for current chromosome
        group_B = groups_B.get_group(chrom)
        val_count = len(group_A)
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