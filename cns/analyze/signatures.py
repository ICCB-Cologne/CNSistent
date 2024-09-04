from cns.utils.assemblies import hg19
import numpy as np
import pandas as pd

# count segments per chromosome and subtract 1
def calc_breaks_per_chr(cns):
    breaks = cns.reset_index().groupby(["sample_id", "chrom"]).size()
    breaks = breaks.reset_index(name="breaks")
    breaks["breaks"] = breaks["breaks"] - 1
    return breaks


def add_breaks_per_sample(cns, samples, assembly=hg19):
    res = samples.copy()
    breaks_per_chr = calc_breaks_per_chr(cns)
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        column = f"breaks_{suffix}"
        res[column] = breaks_per_chr.query("chrom in @names").groupby("sample_id")["breaks"].sum()
        res[column] = res[f"breaks_{suffix}"].fillna(0).astype(int)

    res["breaks_tot"] = res["breaks_aut"] + res["breaks_sex"]
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
        column = f"step_{suffix}"
        res[column] = step_per_chr.query("chrom in @names").groupby("sample_id")["step"].mean()
        res[column] = res[f"step_{suffix}"].fillna(0).astype(float)
    
    res["step_tot"] = (res["step_aut"] * len(assembly.aut_names) + res["step_sex"] * len(assembly.sex_names)) / len(assembly.chr_names)
    return res