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


def calc_breaks_per_sample(cns_df, samples, in_column, assembly=hg19):
    res = samples.copy()
    cns_subset_df = cns_df[["sample_id", "chrom", "start", "end", in_column]]
    merged_cns_df = merge_neighbours(cns_subset_df, in_column)
    breaks_per_chr = calc_breaks_per_chr(merged_cns_df)
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        res[f"breaks_{in_column}_{suffix}"] = (
            breaks_per_chr.query("chrom in @names")
            .groupby("sample_id")["breaks"]
            .sum()
            .reindex(res.index)
            .fillna(0)
            .astype(np.int64)
        )

    res[f"breaks_{in_column}_tot"] = res[f"breaks_{in_column}_aut"] + res[f"breaks_{in_column}_sex"]
    return res


def calc_step_per_chr(cns, col):
    groups = cns.groupby(["sample_id", "chrom"])
    res = []
    for group in groups:
        print(group[1])
        vals = group[1][col].values
        val = 0 if len(vals) < 2 else np.abs(np.diff(vals)).mean()
        res.append((group[0][0], group[0][1], val))
    return pd.DataFrame(res, columns=["sample_id", "chrom", f"step"])


def step_per_sample(cns, samples, assembly=hg19):
    res = samples.copy()
    step_per_chr = calc_step_per_chr(cns, "major_cn")
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        res[f"step_{suffix}"] = (
            step_per_chr.query("chrom in @names")
            .groupby("sample_id")["step"]
            .mean()
            .reindex(res.index)
            .fillna(0)
            .astype(np.int64)
        )

    res["step_tot"] = res["step_aut"] * assembly.aut_count + res["step_sex"] * assembly.sex_count / assembly.chr_count
    return res
