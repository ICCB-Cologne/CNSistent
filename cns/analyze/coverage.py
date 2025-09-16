import numpy as np
import pandas as pd
from cns.utils.conversions import calc_lengths
from cns.utils.selection import only_aut, only_sex
from cns.utils.assemblies import hg19


def get_covered_bases(nan_bases_df, samples_df, either_allele):
    """
    Calculates the number of covered bases for each sample.

    Parameters
    ----------
    nan_bases_df : pandas.DataFrame
        DataFrame containing CNS data with NaN values indicating uncovered bases.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    either_allele : bool
        If True, considers either allele for coverage. If False, considers both alleles.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the number of covered bases for each sample.
    """
    res_df = samples_df.copy()
    label = "any" if either_allele else "both"
    aut_df = only_aut(nan_bases_df)
    aut_df_len = calc_lengths(aut_df)
    sex_df = only_sex(nan_bases_df)
    sex_df_len = calc_lengths(sex_df)
    # Group the differences by sample_id and compute the sum for each group
    res_df[f"cover_{label}_aut"] = (
        aut_df_len.groupby(aut_df["sample_id"]).sum().reindex(res_df.index).fillna(0).astype(np.int64)
    )
    res_df[f"cover_{label}_sex"] = (
        sex_df_len.groupby(sex_df["sample_id"]).sum().reindex(res_df.index).fillna(0).astype(np.int64)
    )
    res_df[f"cover_{label}_all"] = res_df[f"cover_{label}_aut"] + res_df[f"cover_{label}_sex"]
    return res_df


def get_missing_chroms(cns_df, samples_df, segs=None, assembly=hg19):
    """
    Identifies missing chromosomes for each sample.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    segs : dict, optional
        Dictionary of segments to consider. If None, all segments are considered.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the count of chromosomes and missing chromosomes for each sample.
    """
    res = samples_df.copy()
    # create a serise where the value is sex_xy if expected_chrs == 'xy' lese it is sex_xx
    xy_names = assembly.aut_names + ["chrX", "chrY"]
    xx_names = assembly.aut_names + ["chrX"]
    if segs is not None:
        seg_chrs = segs.keys()
        xy_names = [x for x in xy_names if x in seg_chrs]
        xx_names = [x for x in xx_names if x in seg_chrs]

    expected_chrs = res["sex"].map({"xy": xy_names, "xx": xx_names, "NA": xx_names})
    tot_chrs = cns_df.groupby("sample_id")["chrom"].unique()

    merged = pd.DataFrame([expected_chrs, tot_chrs]).T
    diff = merged.apply(lambda x: np.setdiff1d(x.iloc[0], x.iloc[1]), axis=1)

    res["chrom_count"] = tot_chrs.apply(lambda x: len(x)).reindex(res.index).fillna(0).astype(np.int64)
    res["chrom_missing"] = diff
    return res
