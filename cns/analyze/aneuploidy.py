import pandas as pd
from cns.utils.assemblies import hg19
from cns.utils.files import find_cn_cols_if_none


def get_expected_ploidy(column, chrom, is_xy, assembly=hg19):
    if chrom == "chrX":
        if is_xy:
            if column == "minor_cn" or column == "hapY_cn":
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
            if column == "minor_cn" or column == "hapY_cn":
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


# per chromosome and and allele find the number of aneuploid bases
def calc_ane_per_chrom(cns_df, samples_df, cn_column, assembly=hg19):
    aneuploidy = []
    ane_col = f"ane_{cn_column}"
    cns_indexed = cns_df.set_index(["sample_id", "chrom"])
    for (sample, chrom), group in cns_indexed.groupby(level=[0, 1]):
        is_xy = samples_df.loc[sample]["sex"] == "xy"
        sample_data = [sample, chrom]
        values = group[["start", "end", cn_column]].values
        expected_ploidy = get_expected_ploidy(cn_column, chrom, is_xy, assembly)
        aneuploid = values[values[:, 2] != expected_ploidy]
        anu_len = (aneuploid[:,1] - aneuploid[:,0]).sum()
        sample_data.append(anu_len)
        aneuploidy.append(sample_data)
    res = pd.DataFrame(aneuploidy, columns=["sample_id", "chrom"] + [ane_col]).set_index(["sample_id", "chrom"])
    return res


def _calc_ane_per_sample_int(cns_df, samples_df, cn_columns, chroms_list, assembly=hg19):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    sel_df = cns_df.query("chrom in @chroms_list")
    ane_res = [calc_ane_per_chrom(sel_df, samples_df, cn_column, assembly)  for cn_column in cn_columns]
    chr_ane_df = pd.concat(ane_res, axis=1)
    ids = cns_df["sample_id"].unique()
    sample_sum = (
        chr_ane_df
        .groupby("sample_id")
        .apply(lambda g: g.sum())
    )
    # add additional rows to match the ids
    sample_sum = sample_sum.reindex(ids).fillna(0).astype(int)
    return sample_sum


def calc_aut_aneuploidy(cns_df, samples_df, cn_columns=None, assembly=hg19):
    return _calc_ane_per_sample_int(cns_df, samples_df, cn_columns, assembly.aut_names, assembly)


def calc_sex_aneuploidy(cns_df, samples_df, cn_columns=None, assembly=hg19):
    return _calc_ane_per_sample_int(cns_df, samples_df, cn_columns, assembly.sex_names, assembly)


def norm_aut_aneuploidy(samples_df, cn_columns=None, assembly=hg19):
    cn_columns = find_cn_cols_if_none(samples_df, cn_columns)
    res = samples_df.copy()
    for cn_column in cn_columns:        
        res[cn_column + "_frac"] = res.apply(lambda x: x[cn_column] / assembly.aut_len, axis=1)
    return res

def norm_sex_aneuploidy(samples_sex, samples_df, cn_columns=None, assembly=hg19):
    cn_columns = find_cn_cols_if_none(samples_df, cn_columns)
    sex_info = samples_sex[["sex"]]
    merged_df = samples_df.merge(sex_info, left_index=True, right_index=True, how="left")
    
    xx_len = assembly.chr_lens["chrX"]
    xy_len = assembly.chr_lens["chrY"] + xx_len
    
    for cn_column in cn_columns:        
        merged_df["expected_length"] = merged_df.apply(lambda x: xy_len if x["sex"] == "xy" else xx_len, axis=1)
        merged_df[cn_column + "_frac"] = merged_df[cn_column] / merged_df["expected_length"]

    merged_df = merged_df.drop(columns=["sex", "expected_length"])
    return merged_df

def norm_gen_aneuploidy(samples_sex, samples_df, cn_columns=None, assembly=hg19):
    cn_columns = find_cn_cols_if_none(samples_df, cn_columns)
    sex_info = samples_sex[["sex"]]
    merged_df = samples_df.merge(sex_info, left_index=True, right_index=True, how="left")
    xx_len = assembly.aut_len + assembly.chr_lens["chrX"]
    xy_len = assembly.gen_len
    
    for cn_column in cn_columns:        
        merged_df["expected_length"] = merged_df.apply(lambda x: xy_len if x["sex"] == "xy" else xx_len, axis=1)
        merged_df[cn_column + "_frac"] = merged_df[cn_column] / merged_df["expected_length"]

    merged_df = merged_df.drop(columns=["sex", "expected_length"])
    return merged_df

def norm_gen_aneuploidy():
    raise NotImplementedError("Not implemented yet")