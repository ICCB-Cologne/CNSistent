from os.path import join as exists
import pandas as pd

from cns.process.segments import regions_to_segments
from cns.process.imputation import fill_sex_if_missing


def load_cns(path, sort=False, change_coords=True):
    cns = pd.read_csv(path, sep="\t")
    cns = rename_columns(cns)
    cns["start"] = cns["start"].astype(int)
    cns["end"] = cns["end"].astype(int)
    if change_coords:
        cns.loc[:, "start"] -= 1
    if sort:
        cns.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return cns


def save_cns(cns, path, sort=False, change_coords=True):
    cns = cns.copy()
    if change_coords:
        cns.loc[:, "start"] += 1
    if sort:
        cns.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    cns.to_csv(path, sep="\t", index=False)


def load_samples(path):
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    samples_df = pd.read_csv(path, sep="\t", index_col=0)
    samples_df.index.name = "sample_id"
    if "sex" not in samples_df.columns:
        samples_df["sex"] = "NA"     
    return samples_df   


def samples_df_from_cns_df(cns_df, fill_sex=True):
    ids = cns_df["sample_id"].unique()
    samples_df = pd.DataFrame({"sample_id": ids})
    samples_df["sex"] = "NA"
    samples_df.set_index("sample_id", inplace=True)
    if fill_sex:
        samples_df = fill_sex_if_missing(cns_df, samples_df)
    return samples_df


def save_regions(seg_df, path, change_coords=True):
    if change_coords:
        seg_df = seg_df.copy()
        seg_df.loc[:, "start"] += 1
    sel = seg_df[["chrom", "start", "end"]]
    sel.to_csv(path, sep="\t", index=False)


def load_regions(path, change_coords=True):
    if path == "" or path is None:
        return None
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    regs = pd.read_csv(path, sep="\t")
    # check that columns "chrom", "start" and "end" exist, more colums may be present
    if not all([col in regs.columns for col in ["chrom", "start", "end"]]):
        raise ValueError(f"File {path} must have columns 'chrom', 'start' and 'end'.")
    segs = regions_to_segments(regs, change_coords)
    return segs
    

def rename_columns(cns_df):
    res_df = cns_df.copy()
    if res_df.columns.size < 6:
        raise ValueError(
            "cns_df must have first 6 columns in the following order: ",
            "sample_id",
            "chrom",
            "start",
            "end",
            "major_cn",
            "minor_cn",
        )
    # drop columns aftr the 6th
    res_df = res_df.iloc[:, :6]
    res_df.columns = ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"]
    # check if any value in column chrom starts with "chr"
    res_df["chrom"] = res_df["chrom"].astype(str)
    if not res_df.chrom.str.startswith("chr").any():
        res_df.chrom = "chr" + res_df.chrom.astype(str)
    res_df.sort_values(
        by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True
    )
    return res_df


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


def test():
    print("test")