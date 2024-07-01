from os.path import join as exists
import pandas as pd

from cns.utils.conversions import cns_to_segments
from cns.utils.assemblies import hg19


def load_cns(path, cn_col_no=0, sort=False, change_coords=True, no_sample=False, header=True):
    cns_df = pd.read_csv(path, sep="\t", header=0 if header else None)
    if no_sample:
        cns_df["sample"] = "dummy"
    cns_df = canonize_cns_df(cns_df, cn_col_no)
    if change_coords:
        cns_df.loc[:, "start"] -= 1
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return cns_df


def save_cns(cns_df, path, sort=False, change_coords=True, no_sample=False, header=False, write_mode="w"):
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    if change_coords:
        cns_df.loc[:, "start"] += 1

    to_save = cns_df.drop(columns="sample_id") if no_sample else cns_df
    to_save.to_csv(path, sep="\t", index=False, header=header, mode=write_mode)

    if change_coords:
        cns_df.loc[:, "start"] -= 1


def load_samples(path):
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    samples_df = pd.read_csv(path, sep="\t", index_col=0)
    samples_df.index.name = "sample_id"
    if "sex" not in samples_df.columns:
        samples_df["sex"] = "NA"     
    return samples_df   


def fill_sex_if_missing(cns, samples):
    res = samples.copy()
    # Set found_sex to True for each sample if there is chrY, otherwise set it to False
    found_sex = cns.groupby("sample_id").apply(lambda x: "chrY" in x["chrom"].values)
    found_sex = found_sex.map({True: "xy", False: "xx"})
    # replace values in samples["sex"] with found_sex if samples["sex"] is not xy or xx
    res.loc[~res["sex"].isin(["xy", "xx"]), "sex"] = found_sex
    return res


def samples_df_from_cns_df(cns_df, fill_sex=True):
    ids = cns_df["sample_id"].unique()
    samples_df = pd.DataFrame({"sample_id": ids})
    samples_df["sex"] = "NA"
    samples_df.set_index("sample_id", inplace=True)
    if fill_sex:
        samples_df = fill_sex_if_missing(cns_df, samples_df)
    return samples_df


def save_regions(seg_df, path, change_coords=True, header=True, write_mode="w"):
    if change_coords:
        seg_df = seg_df.copy()
        seg_df.loc[:, "start"] += 1
    sel = seg_df[["chrom", "start", "end"]]
    sel.to_csv(path, sep="\t", index=False, header=header, mode=write_mode)


def load_regions(path, change_coords=True, header=True):
    if path == "" or path is None:
        return None
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    regs = pd.read_csv(path, sep="\t", header=0 if header else None)
    # check that columns "chrom", "start" and "end" exist, more colums may be present
    if header:
        if not all([col in regs.columns for col in ["chrom", "start", "end"]]):
            raise ValueError(f"File {path} must have columns 'chrom', 'start' and 'end'.")
    else:
        if len(regs.columns) < 3:
            raise ValueError(f"File {path} must have at least 3 columns.")
        regs.columns = ["chrom", "start", "end"]    
    segs = cns_to_segments(regs, change_coords)
    return segs
    

def canonize_cns_df(cns_df, cn_columns_no=0, assembly=hg19, print_info=False):
    if cn_columns_no > 0:
        if cns_df.columns.size < 4 + cn_columns_no:
            raise ValueError(f"Not enough columns in the CNS file, expected at least 4 + {cn_columns_no}, got {cns_df.columns.size}.")
    else:
        if cns_df.columns.size < 5:
            raise ValueError("Not enough columns in the CNS file, expected at least 5, got {cns_df.columns.size}.")

    # if the column sample_id does not exist, rename the first column to sample_id
    if "sample_id" not in cns_df.columns:
        cns_df.columns = ["sample_id"] + cns_df.columns[1:].tolist()

    # if the column chrom does not exist, rename the second column to chrom
    if "chrom" not in cns_df.columns:
        cns_df.columns = cns_df.columns[:1].tolist() + ["chrom"] + cns_df.columns[2:].tolist()
    chrom_vals = cns_df["chrom"].unique()
    if not any([chrom in assembly.chr_names for chrom in chrom_vals]):
        raise ValueError(f"No chrom found. Chromosome values must be in {assembly.chr_names}, got {chrom_vals}.")
    not_known = [chrom not in assembly.chr_names for chrom in chrom_vals]
    if len(not_known) > 0 and print_info:
        print(f"Found chromosomes not in assembly: {chrom_vals}, these will be dropped.")
    cns_df = cns_df[~cns_df["chrom"].isin(chrom_vals[not_known])]

    # if the column start does not exist, rename the third column to start
    if "start" not in cns_df.columns:
        cns_df.columns = cns_df.columns[:2].tolist() + ["start"] + cns_df.columns[3:].tolist()    
    cns_df["start"] = cns_df["start"].astype(int)

    # if the column end does not exist, rename the fourth column to end
    if "end" not in cns_df.columns:
        cns_df.columns = cns_df.columns[:3].tolist() + ["end"] + cns_df.columns[4:].tolist()
    cns_df["end"] = cns_df["end"].astype(int)
    
    if cn_columns_no <= 0:
        cn_columns = get_cn_columns(cns_df)
        # if cn_columns is empty, take the 5 + columns
        if len(cn_columns) == 0:
            cn_columns = cns_df.columns[5:].tolist()
    else:
        cn_columns = cns_df.columns[4:4 + cn_columns_no].tolist()
        cn_columns_map = {col : (col if is_cn_column(col) else f"{col}_cn") for col in cn_columns}
        cns_df.rename(columns=cn_columns_map, inplace=True)
        cn_columns = list(cn_columns_map.values())
    if print_info:
        print(f"Using CN columns: {cn_columns}")

    sel_columns = ["sample_id", "chrom", "start", "end"] + cn_columns
    cns_df = cns_df[sel_columns]
    return cns_df


def is_cn_column(column):
    if not isinstance(column, str):
        return False
    return column.endswith("cn") or column.endswith("CN") or column.startswith("cn") or column.startswith("CN")


def get_cn_columns(df):
    return [col for col in df.columns if is_cn_column(col)]


def find_cn_cols_if_none(cns_df, cn_cols):
    if cn_cols is None:
        cn_cols = get_cn_columns(cns_df)
    if len(cn_cols) == 0:
        raise ValueError("No CN columns found.")
    return cn_cols


def get_ane_cols_if_none(cns_df, ane_cols=None):
    if ane_cols is None:
        return [col for col in cns_df.columns if col.startswith("ane_") and not col.endswith("_frac")]
    return ane_cols


def dataframe_array_split(df, n_splits):
    """
    Splits a DataFrame into n_splits parts as equally as possible.
    
    Parameters:
    - df: The pandas DataFrame to split.
    - n_splits: The number of parts to split the DataFrame into.
    
    Returns:
    - A list of pandas DataFrame objects.
    """
    # Ensure n_splits is a positive integer
    n_splits = max(int(n_splits), 1)
    
    # Calculate the number of rows in each split
    total_rows = len(df)
    rows_per_split = total_rows // n_splits
    remainder = total_rows % n_splits
    
    # Initialize variables to keep track of the current row and the list of splits
    current_row = 0
    splits = []
    
    for i in range(n_splits):
        # Calculate the number of rows for this split
        # Add one to some of the splits to distribute the remainder
        rows_in_split = rows_per_split + (1 if i < remainder else 0)
        
        # Slice the DataFrame for this split and append to the list
        split_df = df.iloc[current_row:current_row + rows_in_split]
        splits.append(split_df)
        
        # Update the current row index
        current_row += rows_in_split
    
    return splits

