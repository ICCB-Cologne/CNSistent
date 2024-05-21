import pandas as pd
from cns.utils.assemblies import hg19


def get_expected_ploidy(column, chrom, is_xy):
    if chrom == "chrX":
        if is_xy:
            if column == "minor_cn":
                return 0
            else:
                return 1
        else:
            if column == "total_cn":
                return 2
            else:
                return 1
    elif chrom == "chrY":
        if is_xy:
            if column == "minor_cn":
                return 0
            else:
                return 1
        else:
            return 0
    else:
        if column == "total_cn":
            return 2
        else:
            return 1


def _sum_chrom_different(values, expected_ploidy=1):
    return values[values[:, 1] != expected_ploidy, 0].sum()


# per chromosome and and allele find the number of aneuploid bases
def calc_ane_per_chrom(cns, samples_indexed):
    aneuploidy = []
    cns_indexed = cns.set_index(["sample_id", "chrom"])
    for (sample, chrom), group in cns_indexed.groupby(level=[0, 1]):
        is_xy = samples_indexed.loc[sample]["sex"] == "xy"
        sample_data = [sample, chrom]
        for column in ["major_cn", "minor_cn", "total_cn"]:
            values = group[["length", column]].values
            expected_ploidy = get_expected_ploidy(column, chrom, is_xy)
            anu_len = _sum_chrom_different(values, expected_ploidy)
            sample_data.append(anu_len)
        aneuploidy.append(sample_data)
    res = pd.DataFrame(
        aneuploidy, columns=["sample_id", "chrom", "ane_major_cn", "ane_minor_cn", "ane_total_cn"]
    )
    return res


def norm_chrom_aneuploidy(aneuploidies, assembly=hg19):
    res = aneuploidies.copy()
    for column in ["ane_major_cn", "ane_minor_cn", "ane_total_cn"]:
        def chrom_norm(row):
            chrom = row["chrom"]
            return row[column] / assembly.chr_lens[chrom]
        res[column + "_frac"] = res.apply(chrom_norm, axis=1)
    return res


def _sum_cn_columns(group):
    return pd.Series(
        {
            "ane_major_cn": group["ane_major_cn"].sum(),
            "ane_minor_cn": group["ane_minor_cn"].sum(),
            "ane_total_cn": group["ane_total_cn"].sum(),
        }
    )

# TODO: Should add a total fraction
def calc_ane_per_sample(aneuploidies, assembly=hg19):
    ids = aneuploidies.index.unique()

    res = []
    for filter in [assembly.aut_names, assembly.sex_names]:
        selection = aneuploidies["chrom"].isin(filter)
        sample_sum = (
            aneuploidies[selection]
            .groupby("sample_id")
            .apply(_sum_cn_columns)
        )
        # add additional rows to match the ids
        sample_sum = sample_sum.reindex(ids).fillna(0).astype(int)
        res.append(sample_sum)

    return res


def norm_aut_aneuploidy(autosomes_sum, assembly=hg19):
    res = autosomes_sum.copy()
    for column in ['ane_major_cn', 'ane_minor_cn', 'ane_total_cn']:
        res[column + "_frac"] = res.apply(lambda x: x[column] / assembly.aut_len, axis=1)
    return res


def norm_sex_aneuploidy(samples_indexed, sex_chromo_sum, assembly=hg19):
    sex_info = samples_indexed[["sex"]]
    merged_df = sex_chromo_sum.merge(sex_info, left_index=True, right_index=True, how="left")
    xx_len = assembly.chr_lens["chrX"]
    xy_len = assembly.chr_lens["chrY"] + xx_len
    # map xy_len to merged_df
    merged_df["expected_length"] = merged_df.apply(lambda x: xy_len if x["sex"] == "xy" else xx_len, axis=1)
    for column in ["ane_major_cn", "ane_minor_cn", "ane_total_cn"]:
        merged_df[column + "_frac"] = merged_df[column] / merged_df["expected_length"]
    # drop the sex column
    merged_df = merged_df.drop(columns=["sex", "expected_length"])
    return merged_df