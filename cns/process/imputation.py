import numpy as np
import pandas as pd
from cns.utils.files import find_cn_cols_if_none


def get_nan_segs(cns_df):
    nans = cns_df[cns_df.isna().any(axis=1)].copy()
    return nans


def add_tails(cns_df, chr_lengths, print_info=True):
    grouped = cns_df.groupby(["sample_id", "chrom"]).agg({"start": "min", "end": "max"})
    grouped = grouped.rename(columns={"start": "min_start", "end": "max_end"})
    grouped = grouped.reset_index()
    missing_ranges = []
    for _, row in grouped.iterrows():
        if row.min_start > 0:
            missing_ranges.append(
                {
                    "sample_id": row.sample_id,
                    "chrom": row.chrom,
                    "start": 0,
                    "end": row.min_start
                }
            )
        if row.max_end < chr_lengths[f"{row.chrom}"]:
            missing_ranges.append(
                    {
                        "sample_id": row.sample_id,
                        "chrom":  row.chrom,
                        "start": row.max_end,
                        "end": chr_lengths[str(row.chrom)]
                    }
            )

    if len(missing_ranges) == 0:
        if print_info:
            print(f"No missing ends found.")
        return cns_df.copy()
    else:
        if print_info:
            print(f"Adding {len(missing_ranges)} missing ends")
        new_cns_df = pd.DataFrame(missing_ranges, columns=cns_df.columns)
        res_df = pd.concat([cns_df, new_cns_df])
        res_df.sort_values(
            by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True
        )
        return res_df


def fill_gaps(cns_df, print_info=True):
    # Iterate over the rows
    new_rows = []
    for i in range(len(cns_df) - 1):
        # Check if the next row has the same 'sample_id' and 'chrom'
        if (
            cns_df.at[i, "sample_id"] == cns_df.at[i + 1, "sample_id"]
            and cns_df.at[i, "chrom"] == cns_df.at[i + 1, "chrom"]
        ):
            # Calculate the range
            range_start = cns_df.at[i, "end"]
            range_end = cns_df.at[i + 1, "start"]
            # If the range is greater than 1, add a new row
            if range_end > range_start:
                new_rows.append(
                    {
                        "sample_id": cns_df.at[i, "sample_id"],
                        "chrom": cns_df.at[i, "chrom"],
                        "start": range_start,
                        "end": range_end,
                    }
                )

    if len(new_rows) == 0:
        if print_info:
            print(f"No gaps found.")
        return cns_df.copy()
    else:
        # Concatenate the cns_dfs
        if print_info:
            print(f"Filling {len(new_rows)} gaps.")
        new_cns_df = pd.DataFrame(new_rows, columns=cns_df.columns)
        res_df = pd.concat([cns_df, new_cns_df])
        res_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
        return res_df


# Add fully missing chromosomes
def add_missing(cns_df, samples_df, chr_lengths, print_info=True):
    res_df = cns_df.set_index("sample_id")
    chromosomes = chr_lengths.keys()

    new_entries = []
    for sample in res_df.index.unique():
        cns_sample_df = res_df.loc[sample]
        sample_chroms = cns_sample_df["chrom"].values
        for chromosome in chromosomes:
            if chromosome not in sample_chroms and (chromosome != "chrY" or samples_df.loc[sample].sex == "xy"):
                new_entry = {
                    "sample_id": sample,
                    "chrom": chromosome,
                    "start": 0,
                    "end": chr_lengths[chromosome],
                }
                new_entries.append(new_entry)

    if len(new_entries) == 0:
        if print_info:
            print(f"No missing chromosomes found.")
        return cns_df.copy()
    else:
        if print_info:
            print(f"Adding {len(new_entries)} missing chromosomes.")
        empty_chrs_df = pd.DataFrame(new_entries)
        res_df.reset_index(inplace=True)
        res_df = pd.concat([res_df, empty_chrs_df])
        res_df.sort_values(
            by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True
        )
        return res_df


# Makes sure that the columns are of the correct type
def _are_mergeable(a, b, cn_columns):
    return (
        a.sample_id == b.sample_id
        and a.chrom == b.chrom
        and a.end == b.start
        and all([(a[col] == b[col]) or (np.isnan(a[col]) and np.isnan(b[col])) for col in cn_columns])
    )


def merge_neighbours(cns_df, cn_columns=None, print_info=True):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)    
    res_df = cns_df.copy()
    idx_to_remove = []

    for i, (index, row) in enumerate(res_df.iterrows()):
        if i != 0 and _are_mergeable(prev, row, cn_columns):
            idx_to_remove.append(i - 1)
            res_df.at[index, "start"] = prev.start
            row.start = prev.start  # update the comparison copy too
        prev = row

    if print_info:
        print(f"Merged entries: {len(idx_to_remove)}")

    # remove from cns_df where idx_to_remove is in the index
    res_df = res_df.drop(res_df.index[idx_to_remove]).sort_values(by=["sample_id", "chrom", "start"], ignore_index=True)
    return res_df


def _is_same_chrom(df, i, j):
    return df.at[j, "sample_id"] == df.at[i, "sample_id"] and df.at[j, "chrom"] == df.at[i, "chrom"]


def create_imputed_entries(cns_df, cn_columns=None, print_info=True):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)    
    new_entries = []
    for i in range(len(cns_df)):
        if np.isnan(cns_df.at[i, cn_columns[0]]) or np.isnan(cns_df.at[i, cn_columns[1]]):
            prev_vals = [cns_df.at[i, cn_columns[k]] for k in range(len(cn_columns))]
            next_vals = list(prev_vals)
            for k in range(len(cn_columns)):
                if np.isnan(prev_vals[k]):
                    j = i - 1
                    while j >= 0 and np.isnan(prev_vals[k]) and _is_same_chrom(cns_df, i, j):
                        prev_vals[k] = cns_df.at[j, cn_columns[k]]
                        j -= 1
                if np.isnan(next_vals[k]):
                    j = i + 1
                    while j < len(cns_df) and np.isnan(next_vals[k]) and _is_same_chrom(cns_df, i, j):
                        next_vals[k] = cns_df.at[j, cn_columns[k]]
                        j += 1
            id = cns_df.at[i, "sample_id"]
            chrom = cns_df.at[i, "chrom"]            
            start = cns_df.at[i, "start"]
            end = cns_df.at[i, "end"]
            if i == 0 or not _is_same_chrom(cns_df, i, i-1):
                new_start = [id, chrom, start, end] + next_vals
                new_entries.append(new_start)
            elif i == len(cns_df) - 1 or not _is_same_chrom(cns_df, i, i+1):
                new_end = [id, chrom, start, end] + prev_vals
                new_entries.append(new_end)
            elif all([prev_vals[k] == next_vals[k] for k in range(len(cn_columns))]):
                new_simple = [id, chrom, start, end] + prev_vals
                new_entries.append(new_simple)
            else:
                midpoint = start + (end - start) // 2
                first_half = [id, chrom, start, midpoint] + prev_vals
                second_half = [id, chrom, midpoint, end] + next_vals
                new_entries.append(first_half)                    
                new_entries.append(second_half)                
        
    imputataion_df = pd.DataFrame(new_entries, columns=cns_df.columns)
    query = ' or '.join([f"{col}.isnull()" for col in cn_columns])
    idx_to_remove = cns_df.query(query).index

    if print_info:
        print(f"New entries: {imputataion_df.shape[0]}")
        print(f"Removed entries: {len(idx_to_remove)}")
    # remove from cns_df where idx_to_remove is in the index
    filtered_df = cns_df.drop(idx_to_remove)
    # concat the new_table to cns_df
    res_df = pd.concat([filtered_df, imputataion_df])
    # sort cns_df by sample_id, chr, start
    res_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return res_df


def fill_nans_with_zeros(cns_df, cn_columns=None, print_info=True):    
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)    
    res_df = cns_df.copy()
    if print_info:
        print(f"Filling {res_df[cn_columns].isna().any(axis=1).sum()} NaN rows with zero")
    # Fully missing chromosomes filled with 0
    for col in cn_columns:
        res_df[col] = res_df[col].fillna(0).astype(int)
    return res_df