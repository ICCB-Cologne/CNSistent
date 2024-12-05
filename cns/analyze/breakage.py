from cns.utils.selection import get_chr_sets
from cns.utils.assemblies import hg19
import numpy as np
import pandas as pd


# count segments per chromosome and subtract 1
def calc_breaks_per_chr(cns_df):
    def con_match_count(group):
        shifted_group = group.shift(-1)
        return (group['end'] == shifted_group['start']).sum()
    return cns_df.groupby(['sample_id', 'chrom']).apply(con_match_count).reset_index(name='breaks')


def calc_breaks_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
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
    def sum_abs_diff(group):
        shifted_group = group.shift(-1)
        consecutive = group['end'] == shifted_group['start']
        abs_diff = (group[cn_col] - shifted_group[cn_col]).abs()
        step_sum = abs_diff[consecutive].sum()
        step_count = consecutive.sum()
        return pd.Series({'step': step_sum, 'count': step_count})
    
    res_df = cns_df.groupby(['sample_id', 'chrom']).apply(sum_abs_diff).reset_index()
    return res_df


def calc_step_per_sample(cns_df, samples_df, cn_col, assembly=hg19):
    res = samples_df.copy()
    step_per_chr = calc_step_per_chr(cns_df, cn_col)
    chrom_types = get_chr_sets(cns_df, assembly)

    for suffix, names in chrom_types.items():
        res[f"step_{cn_col}_{suffix}"] = (
            step_per_chr.query("chrom in @names").groupby("sample_id")
            .apply(lambda x: x["step"].sum() / x["count"].sum() if x["count"].sum() != 0 else 0)
            .reset_index(name="cnstep").set_index("sample_id")
        )
    return res.fillna(0)