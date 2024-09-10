import numpy as np
import pandas as pd
from cns.process.imputation import merge_neighbours
from cns.utils.assemblies import hg19


def calc_step_per_chr(cns_df, cn_col):
    groups = cns_df.groupby(["sample_id", "chrom"])
    res = []
    for group in groups:
        vals = group[1][cn_col].values
        val = 0 if len(vals) < 2 else np.abs(np.diff(vals)).sum()
        res.append((group[0][0], group[0][1], val, len(vals) - 1))
    return pd.DataFrame(res, columns=["sample_id", "chrom", "step", "count"])


def get_step_per_sample(cns_df, samples, cn_col, assembly=hg19):
    res = samples.copy()
    cns_subset_df = cns_df[["sample_id", "chrom", "start", "end", cn_col]]
    merged_cns_df = merge_neighbours(cns_subset_df, cn_col, False)
    step_per_chr = calc_step_per_chr(merged_cns_df, cn_col)
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names, "tot": assembly.chr_names}

    for suffix, names in chrom_types.items():
        res[f"step_{cn_col}_{suffix}"] = (
            step_per_chr.query("chrom in @names").groupby("sample_id")
            .apply(lambda x: x["step"].sum() / x["count"].sum() if x["count"].sum() != 0 else 0)
            .reset_index(name="cnstep").set_index("sample_id")
        )
    return res.fillna(0)
