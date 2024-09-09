import numpy as np
import pandas as pd
from cns.utils.assemblies import hg19


def get_col_changepoint(cns, cn_col):
    groups = cns.groupby("chrom")
    res = []
    for group in groups:
        vals = group[1][cn_col].values
        if len(vals) < 2:
            continue
        res.append(np.abs(np.diff(vals)).mean())
    return np.mean(res)


def calc_step_per_chr(cns, cn_col):
    groups = cns.groupby(["sample_id", "chrom"])
    res = []
    for group in groups:
        vals = group[1][cn_col].values
        val = 0 if len(vals) < 2 else np.abs(np.diff(vals)).mean()
        res.append((group[0][0], group[0][1], val))
    return pd.DataFrame(res, columns=["sample_id", "chrom", f"step"])


def get_step_per_sample(cns, samples, cn_col, assembly=hg19):
    res = samples.copy()
    step_per_chr = calc_step_per_chr(cns, "major_cn")
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        res[f"step_{cn_col}_{suffix}"] = (
            step_per_chr.query("chrom in @names")
            .groupby("sample_id")["step"]
            .sum()
            .reindex(res.index)
            .fillna(0)
            .astype(np.int64)
        )

    res["step_{in_col}_tot"] = (res["step_{in_col}_aut"] + res["step_{in_col}_sex"]) / assembly.chr_count
    res["step_{in_col}_aut"] /= assembly.aut_count
    res["step_{in_col}_sex"] /= assembly.sex_count
    return res
