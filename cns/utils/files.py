from os.path import abspath, exists, join
import pandas as pd

from cns.utils.canonization import canonize_cns_df, is_canonical_cns_df
from cns.utils.conversions import cns_to_segments


def load_cns(path, canonize=False, cn_columns=None, sort=False, change_coords=True, no_sample=False, header=True):
    cns_df = pd.read_csv(path, sep="\t", header=0 if header else None)
    if no_sample:
        cns_df["sample_id"] = "sample_id"
    if canonize or not header:
        cns_df = canonize_cns_df(cns_df, cn_columns)
    elif not is_canonical_cns_df(cns_df):
        raise ValueError("CNS file is not canonical, call load_cns(..., canonize=True, ...) instead.")
    if change_coords:
        cns_df.loc[:, "start"] -= 1
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return cns_df


def save_cns(cns_df, path, sort=False, change_coords=True, no_sample=False, header=True, write_mode="w"):
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


def save_samples(samples_df, path, no_sample=False, header=True):
    if no_sample:
        samples_df = samples_df.reset_index(drop=True)
    samples_df.to_csv(path, sep="\t", index=True, header=header)


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


def save_regions(seg_df, path, change_coords=True, header=True):
    if change_coords:
        seg_df = seg_df.copy()
        seg_df.loc[:, "start"] += 1
    sel = seg_df[["chrom", "start", "end"]]
    sel.to_csv(path, sep="\t", index=False, header=header)


def load_regions(path, change_coords=True, header=True):
    if path == "" or path is None:
        return None
    path = abspath(path)
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
