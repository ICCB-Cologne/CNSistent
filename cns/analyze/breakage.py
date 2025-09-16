from cns.utils.selection import get_chr_sets
from cns.utils.assemblies import hg19
import numpy as np
import pandas as pd


# count segments per chromosome and subtract 1
def calc_breaks_per_chr(cns_df):
    """
    Calculates the number of breakpoints per chromosome for each sample.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the number of breakpoints per chromosome for each sample.
    """
    def con_match_count(group):
        shifted_group = group.shift(-1)
        return (group['end'] == shifted_group['start']).sum()
    # Use groupby object to select only the columns needed, avoiding deprecated behavior
    return cns_df.groupby(['sample_id', 'chrom'])[['end', 'start']].apply(
        lambda group: con_match_count(group.assign(start_next=group['start']))
    ).reset_index(name='breaks')


def calc_breaks_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
    """
    Calculates the number of breakpoints per sample.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    cn_col : str
        Column name for copy number data.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the number of breakpoints per sample.
    """
    res = samples_df.copy()
    breaks_per_chr = calc_breaks_per_chr(cns_df)
    chrom_types = get_chr_sets(cns_df, assembly)

    for suffix, names in chrom_types.items():
        res[f"breaks_{cn_col}_{suffix}"] = (
            breaks_per_chr.query("chrom in @names")
            .groupby("sample_id")["breaks"]
            .sum()
            .reindex(res.index)
            .fillna(0)
            .astype(np.int64)
        )
    return res


def calc_step_per_chr(cns_df, cn_col):
    """
    Calculates the step size per chromosome for each sample.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_col : str
        Column name for copy number data.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the step size per chromosome for each sample.
    """
    def sum_abs_diff(group):
        shifted_group = group.shift(-1)
        consecutive = group['end'] == shifted_group['start']
        abs_diff = (group[cn_col] - shifted_group[cn_col]).abs()
        step_sum = abs_diff[consecutive].sum()
        step_count = consecutive.sum()
        return pd.Series({'step': step_sum, 'count': step_count})
    # Only pass the columns needed to avoid deprecated behavior
    res_df = cns_df.groupby(['sample_id', 'chrom'])[['end', 'start', cn_col]].apply(sum_abs_diff).reset_index()
    return res_df


def calc_step_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
    """
    Calculates the step size per sample.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    cn_col : str
        Column name for copy number data.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the step size per sample.
    """
    res = samples_df.copy()
    step_per_chr = calc_step_per_chr(cns_df, cn_col)
    chrom_types = get_chr_sets(cns_df, assembly)

    for suffix, names in chrom_types.items():
        grouped = step_per_chr.query("chrom in @names").groupby("sample_id")[["step", "count"]]
        cnstep = grouped.sum()
        cnstep["cnstep"] = cnstep.apply(lambda x: x["step"] / x["count"] if x["count"] != 0 else 0, axis=1)
        res[f"step_{cn_col}_{suffix}"] = cnstep["cnstep"]
    return res.fillna(0)