import pandas as pd
from cns.utils.assemblies import hg19
from cns.utils.files import get_ane_cols_if_none, find_cn_cols_if_none


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
def calc_ane_per_chrom(cns_df, samples_indexed, cn_columns=None):
    aneuploidy = []
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    new_columns = [f"ane_{column}" for column in cn_columns]
    cns_indexed = cns_df.set_index(["sample_id", "chrom"])
    for (sample, chrom), group in cns_indexed.groupby(level=[0, 1]):
        is_xy = samples_indexed.loc[sample]["sex"] == "xy"
        sample_data = [sample, chrom]
        for column in cn_columns:
            if column not in group.columns:
                continue
            values = group[["length", column]].values
            expected_ploidy = get_expected_ploidy(column, chrom, is_xy)
            anu_len = _sum_chrom_different(values, expected_ploidy)
            sample_data.append(anu_len)
        aneuploidy.append(sample_data)
    res = pd.DataFrame(aneuploidy, columns=["sample_id", "chrom"] + new_columns)
    return res


def _sum_cn_columns(group, columns):
    res = {}
    for column in columns:
        if column in group.columns:
            res[column] = group[column].sum()
    return pd.Series(res)
    

def calc_ane_per_sample(cns_df, ane_columns=None, assembly=hg19):
    ids = cns_df["sample_id"].unique()
    ane_columns = get_ane_cols_if_none(cns_df, ane_columns)

    res = []
    for filter in [assembly.aut_names, assembly.sex_names]:
        selection = cns_df["chrom"].isin(filter)
        sample_sum = (
            cns_df[selection]
            .groupby("sample_id")
            .apply(lambda g: _sum_cn_columns(g, ane_columns))
        )
        # add additional rows to match the ids
        sample_sum = sample_sum.reindex(ids).fillna(0).astype(int)
        res.append(sample_sum)

    return res


def norm_aut_aneuploidy(autosomes_sum, ane_columns=None, assembly=hg19):
    res = autosomes_sum.copy()
    ane_columns = get_ane_cols_if_none(res, ane_columns)

    for column in ane_columns:
        if column not in res.columns:
            continue
        res[column + "_frac"] = res.apply(lambda x: x[column] / assembly.aut_len, axis=1)
    return res


def norm_sex_aneuploidy(samples_indexed, sex_chrom_sum, ane_columns=None, assembly=hg19):
    ane_columns = get_ane_cols_if_none(sex_chrom_sum, ane_columns)

    sex_info = samples_indexed[["sex"]]
    merged_df = sex_chrom_sum.merge(sex_info, left_index=True, right_index=True, how="left")
    xx_len = assembly.chr_lens["chrX"]
    xy_len = assembly.chr_lens["chrY"] + xx_len
    # map xy_len to merged_df
    merged_df["expected_length"] = merged_df.apply(lambda x: xy_len if x["sex"] == "xy" else xx_len, axis=1)
    for column in ane_columns:
        if column not in merged_df.columns:
            continue
        merged_df[column + "_frac"] = merged_df[column] / merged_df["expected_length"]
    # drop the sex column
    merged_df = merged_df.drop(columns=["sex", "expected_length"])
    return merged_df