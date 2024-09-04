from cns.utils.assemblies import hg19

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

