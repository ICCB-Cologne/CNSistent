from cns.process.imputation import merge_neighbours
from cns.utils.assemblies import hg19
import numpy as np
import pandas as pd


# count segments per chromosome and subtract 1
def calc_breaks_per_chr(cns):
    breaks = cns.reset_index().groupby(["sample_id", "chrom"]).size()
    breaks = breaks.reset_index(name="breaks")
    breaks["breaks"] = breaks["breaks"] - 1
    return breaks


def calc_breaks_per_sample(cns_df, samples, cn_col, assembly=hg19):
    res = samples.copy()
    cns_subset_df = cns_df[["sample_id", "chrom", "start", "end", cn_col]]
    merged_cns_df = merge_neighbours(cns_subset_df, cn_col, False)
    breaks_per_chr = calc_breaks_per_chr(merged_cns_df)
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        res[f"breaks_{cn_col}_{suffix}"] = (
            breaks_per_chr.query("chrom in @names")
            .groupby("sample_id")["breaks"]
            .sum()
            .reindex(res.index)
            .fillna(0)
            .astype(np.int64)
        )

    res[f"breaks_{cn_col}_tot"] = res[f"breaks_{cn_col}_aut"] + res[f"breaks_{cn_col}_sex"]
    return res
