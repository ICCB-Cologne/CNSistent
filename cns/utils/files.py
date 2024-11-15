from os.path import abspath, exists
import numpy as np
import pandas as pd

from cns.utils.canonization import canonize_cns_df, canonize_sample_id, is_canonical_cns_df
from cns.utils.conversions import df_to_segs, segs_to_df
from cns.utils.logging import log_warn
from cns.utils.assemblies import hg19


def load_cns(path, canonize=False, cn_columns=None, sep="\t", sort=False, change_coords=True, assembly=hg19, print_info=False):
    cns_df = pd.read_csv(path, sep=sep, low_memory=False)
    if canonize:
        cns_df = canonize_cns_df(cns_df, cn_columns, False, assembly, print_info)
    elif not is_canonical_cns_df(cns_df):
        raise ValueError("CNS file is not canonical, call load_cns(..., canonize=True, ...) instead.")
    if change_coords:
        cns_df.loc[:, "start"] -= 1
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return cns_df


def save_cns(cns_df, path, sort=False, change_coords=True, mode="w"):
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    if change_coords:
        cns_df.loc[:, "start"] += 1
    cns_df.to_csv(path, sep="\t", index=False, mode=mode, header=mode=="w")
    if change_coords:
        cns_df.loc[:, "start"] -= 1


def load_samples(path, cananonize=False, sep="\t"):
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    samples_df = pd.read_csv(path, sep=sep)
    if cananonize:
        samples_df = canonize_sample_id(samples_df)
    elif "sample_id" not in samples_df.columns:
        raise ValueError("Column 'sample_id' not found in the samples file. Call load_samples(..., cananonize=True, ...) instead.")	
    if "sex" not in samples_df.columns:
        samples_df["sex"] = "NA"     
    samples_df.set_index("sample_id", inplace=True)
    return samples_df   


def save_samples(samples_df, path, mode='w'):
    samples_df.to_csv(path, sep="\t", index=True, mode=mode, header=mode=="w")


def fill_sex_if_missing(cns_df, samples_df):
    res_df = samples_df.copy()
    # Set found_sex to True for each sample if there is chrY, otherwise set it to False
    found_sex = cns_df.groupby("sample_id").apply(lambda x: "chrY" in x["chrom"].values)
    found_sex = found_sex.map({True: "xy", False: "xx"})
    if "sex" in res_df.columns:
        res_df["found_sex"] = found_sex
        condition = (res_df["sex"] == "xx") & (res_df["found_sex"] == "xy")
        indices = res_df[condition].index
        if len(indices) > 0:
            log_warn(f"Found samples where sex is xx in data but chrY has CNs assigned: {indices.tolist()}. "\
                    "This may result in an incorrect proportions of sex-chromosome features.")          
        res_df.drop(columns=["found_sex"], inplace=True)      
    # replace values in samples["sex"] with found_sex if samples["sex"] is not xy or xx
    res_df.loc[~res_df["sex"].isin(["xy", "xx"]), "sex"] = found_sex
    return res_df


def find_y_column(cns_df, samples_df, cn_columns):
    res_df = samples_df.copy()
    res_df["y_col"] = "NA"
    for sample_id, group_df in cns_df.groupby("sample_id"):
        chr_y = group_df.query("chrom == 'chrY'")
        if len(chr_y) > 0:
            if (chr_y[cn_columns[1]] == 0).all():
                res_df.loc[sample_id, "y_col"] = cn_columns[0]
            elif (chr_y[cn_columns[0]] == 0).all():
                res_df.loc[sample_id, "y_col"] = cn_columns[1]
            else:
                log_warn(f"Sample {sample_id} has non-zero chrY CNs for both haplotypes.") 
    return res_df


def samples_df_from_cns_df(cns_df, fill_sex=True):
    ids = cns_df["sample_id"].unique()
    samples_df = pd.DataFrame({"sample_id": ids})
    samples_df["sex"] = "NA"
    samples_df.set_index("sample_id", inplace=True)
    if fill_sex:
        samples_df = fill_sex_if_missing(cns_df, samples_df)
    return samples_df


def save_segments(segs, path):    
    is_bed = path.lower().endswith(".bed")
    if not is_bed:
        log_warn(f"Segments file {path} is not bed file, the coordinates will be 1-based.")
    seg_df = segs_to_df(segs)
    if not is_bed:
        seg_df = seg_df.copy()
        seg_df.loc[:, "start"] += 1
    sel = seg_df[["chrom", "start", "end", "name"]]
    sel.to_csv(path, sep="\t", index=False, header=not is_bed)


def load_segments(path):
    is_bed = path.lower().endswith(".bed")
    if not is_bed:
        log_warn(f"Segments file {path} is not bed file, the coordinates will be 1-based.")
    if path == "" or path is None:
        return None
    path = abspath(path)
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    segs_df = pd.read_csv(path, sep="\t", header=(None if is_bed else 0))
    # check that columns "chrom", "start" and "end" exist, more colums may be present
    if not is_bed:
        if not all([col in segs_df.columns for col in ["chrom", "start", "end"]]):
            raise ValueError(f"File {path} must have columns 'chrom', 'start' and 'end'.")
        if "name" not in segs_df.columns:
            segs_df["name"] = np.arange(len(segs_df))
    else:
        if len(segs_df.columns) < 3:
            raise ValueError(f"File {path} must have at least 3 columns.")
        elif len(segs_df.columns) == 3:            
            segs_df["name"] = np.arange(len(segs_df))
        elif len(segs_df.columns) > 4:
            print(f"Warning: File {path} has more than 4 columns. Only the first 4 columns are used.")
            segs_df = segs_df.iloc[:, :4]                
        segs_df.columns = ["chrom", "start", "end", "name"]  
    if not is_bed:
        segs_df.loc[:, "start"] -= 1
    if len(segs_df.columns) == 3:
        segs_df["name"] = np.arange(len(segs_df))

    return df_to_segs(segs_df)
    

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
        if rows_in_split > 0:
            splits.append(split_df)            
            # Update the current row index
            current_row += rows_in_split
    
    return splits
