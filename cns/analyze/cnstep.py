import numpy as np
import pandas as pd
from cns.utils.assemblies import hg19


def get_col_changepoint(cns, col):
	groups = cns.groupby("chrom")
	res = []
	for group in groups:
		vals = group[1][col].values
		if len(vals) < 2:
			continue
		res.append(np.abs(np.diff(vals)).mean())
	return np.mean(res)


def calc_step_per_chr(cns, col):
    groups = cns.groupby(["sample_id", "chrom"])
    res = []
    for group in groups:
        vals = group[1][col].values
        val = 0 if len(vals) < 2 else np.abs(np.diff(vals)).mean()   
        res.append((group[0][0], group[0][1], val))
    return pd.DataFrame(res, columns=["sample_id", "chrom", f"step"])
    

def step_per_sample(cns, samples, cn_col, assembly=hg19):
    res = samples.copy()
    step_per_chr = calc_step_per_chr(cns, cn_col)
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, _ in chrom_types.items():
        step_col = f"step_{cn_col}_{suffix}"
        res[step_col] = step_per_chr.query("chrom in @names").groupby("sample_id")["step"].mean()
        res[step_col] = res[f"step_{suffix}"].fillna(0).astype(float)
    
    res["step_tot"] = (res["step_aut"] * len(assembly.aut_names) + res["step_sex"] * len(assembly.sex_names)) / len(assembly.chr_names)
    return res


def get_cn_steps_per_column(cns, samples, cn_col, assembly=hg19):
    res = samples.copy()
    step_per_chr = calc_step_per_chr(cns, "major_cn")
    chrom_types = {"aut": assembly.aut_names, "sex": assembly.sex_names}

    for suffix, names in chrom_types.items():
        cn_col = f"step_{suffix}"
        res[cn_col] = step_per_chr.query("chrom in @names").groupby("sample_id")["step"].mean()
        res[cn_col] = res[f"step_{suffix}"].fillna(0).astype(float)
    
    res["step_tot"] = (res["step_aut"] * len(assembly.aut_names) + res["step_sex"] * len(assembly.sex_names)) / len(assembly.chr_names)
    return res