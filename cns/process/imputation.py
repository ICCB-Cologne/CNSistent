import numpy as np
import pandas as pd
from cns.utils.canonization import find_cn_cols_if_none
from cns.utils.logging import log_info


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
        log_info(print_info, f"No missing ends found.")
        return cns_df.copy()
    else:
        log_info(print_info, f"Adding {len(missing_ranges)} missing ends")
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
        log_info(print_info, f"No gaps found.")
        return cns_df.copy()
    else:
        # Concatenate the cns_dfs
        log_info(print_info, f"Filling {len(new_rows)} gaps.")
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
        log_info(print_info, f"No missing chromosomes found.")
        return cns_df.copy()
    else:
        log_info(print_info, f"Adding {len(new_entries)} missing chromosomes.")
        empty_chrs_df = pd.DataFrame(new_entries)
        res_df.reset_index(inplace=True)
        res_df = pd.concat([res_df, empty_chrs_df])
        res_df.sort_values(
            by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True
        )
        return res_df


def remove_outliers(cns_df, chr_lengths, print_info=True):
    res_df = cns_df.copy()
    idx_to_remove = []
    for i, (index, row) in enumerate(res_df.iterrows()):
        if row.start < 0:
            if row.end < 0:
                idx_to_remove.append(i)
                log_info(print_info, f"Removed outlier:\n{row}")
            else:
                res_df.at[index, "start"] = 0
        if row.end > chr_lengths[row.chrom]:
            if row.start > chr_lengths[row.chrom]:
                idx_to_remove.append(i)
                log_info(print_info, f"Removed outlier:\n{row}")
            else:
                res_df.at[index, "end"] = chr_lengths[row.chrom]

    log_info(print_info, f"Removed outliers: {len(idx_to_remove)}")
    # remove from cns_df where idx_to_remove is in the index
    res_df = res_df.drop(res_df.index[idx_to_remove]).sort_values(by=["sample_id", "chrom", "start"], ignore_index=True)
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

    log_info(print_info, f"Merged entries: {len(idx_to_remove)}")

    # remove from cns_df where idx_to_remove is in the index
    res_df = res_df.drop(res_df.index[idx_to_remove]).sort_values(by=["sample_id", "chrom", "start"], ignore_index=True)
    return res_df


def _is_same_chrom(df, i, j):
    return df.at[j, "sample_id"] == df.at[i, "sample_id"] and df.at[j, "chrom"] == df.at[i, "chrom"]


def _impute_extend(cns_df, cn_columns, print_info=True):
    new_entries = []
    for i in range(len(cns_df)):
        if any(np.isnan(cns_df.at[i, col]) for col in cn_columns):
            prev_vals = [np.nan for _ in cn_columns]
            next_vals = list(prev_vals)
            for k in range(len(cn_columns)):
                j = i
                while j >= 0 and np.isnan(prev_vals[k]) and _is_same_chrom(cns_df, i, j):
                    prev_vals[k] = cns_df.at[j, cn_columns[k]]
                    j -= 1
                j = i
                while j < len(cns_df) and np.isnan(next_vals[k]) and _is_same_chrom(cns_df, i, j):
                    next_vals[k] = cns_df.at[j, cn_columns[k]]
                    j += 1
                if np.isnan(prev_vals[k]) and np.isnan(next_vals[k]):
                    prev_vals[k] = 0
                    next_vals[k] = 0
                elif np.isnan(prev_vals[k]):
                    prev_vals[k] = next_vals[k]
                elif np.isnan(next_vals[k]):
                    next_vals[k] = prev_vals[k]
            
            id = cns_df.at[i, "sample_id"]
            chrom = cns_df.at[i, "chrom"]            
            start = cns_df.at[i, "start"]
            end = cns_df.at[i, "end"]
            if all([prev_vals[k] == next_vals[k] for k in range(len(cn_columns))]):
                new_simple = [id, chrom, start, end] + prev_vals
                new_entries.append(new_simple)
            else:
                midpoint = start + (end - start) // 2
                first_half = [id, chrom, start, midpoint] + prev_vals
                second_half = [id, chrom, midpoint, end] + next_vals
                new_entries.append(first_half)                    
                new_entries.append(second_half)                
    
    new_cols = ["sample_id", "chrom", "start", "end"] + cn_columns
    imputation_df = pd.DataFrame(new_entries, columns=new_cols)
    query = ' or '.join([f"{col}.isnull()" for col in cn_columns])
    idx_to_remove = cns_df.query(query).index

    log_info(print_info, f"New entries: {imputation_df.shape[0]}\nRemoved entries: {len(idx_to_remove)}")
    # remove from cns_df where idx_to_remove is in the index
    filtered_df = cns_df.drop(idx_to_remove)
    # concat the new_table to cns_df
    res_df = pd.concat([filtered_df, imputation_df])
    # sort cns_df by sample_id, chr, start
    res_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return res_df


def _impute_diploid(cns_df, samples_df, cn_columns, print_info=True):
    if len(cn_columns) > 2 or len(cn_columns) < 1:
        raise ValueError("Diploid imputation can only be done for one (total CN) or two (major, minor CN) columns.")
    
    aut_df = cns_df.query("chrom != 'chrX' and chrom != 'chrY'")
    xx_samples = samples_df.query("sex == 'xx'").index
    xy_samples = samples_df.query("sex == 'xy'").index
    xx_cns_df = cns_df[cns_df["sample_id"].isin(xx_samples)]
    xx_x_chrom_df = xx_cns_df.query("chrom == 'chrX'")
    xx_y_chrom_df = xx_cns_df.query("chrom == 'chrY'")
    xy_cns_df = cns_df[cns_df["sample_id"].isin(xy_samples)]
    xy_x_chrom_df = xy_cns_df.query("chrom == 'chrX'")
    xy_y_chrom_df = xy_cns_df.query("chrom == 'chrY'")
    if len(cn_columns) == 2:
        aut_df[cn_columns] = aut_df[cn_columns].fillna(1)
        xx_x_chrom_df[cn_columns] = xx_x_chrom_df[cn_columns].fillna(1)      
        xx_y_chrom_df[cn_columns] = xx_y_chrom_df[cn_columns].fillna(0)
        xy_x_chrom_df[cn_columns[0]] = xy_x_chrom_df[cn_columns[0]].fillna(1)
        xy_x_chrom_df[cn_columns[1]] = xy_x_chrom_df[cn_columns[1]].fillna(0)
        xy_y_chrom_df[cn_columns[0]] = xy_y_chrom_df[cn_columns[0]].fillna(1)
        xy_y_chrom_df[cn_columns[0]] = xy_y_chrom_df[cn_columns[0]].fillna(0)
    else:
        col = cn_columns[0]
        aut_df[col] = aut_df[col].fillna(2)
        xx_x_chrom_df[col] = xx_x_chrom_df[col].fillna(2)
        xx_y_chrom_df[col] = xx_y_chrom_df[col].fillna(0)
        xy_x_chrom_df[col] = xy_x_chrom_df[col].fillna(1)
        xy_y_chrom_df[col] = xy_y_chrom_df[col].fillna(1)

        # Update the original dataframe in place
    cns_df.update(aut_df)
    cns_df.update(xx_x_chrom_df)
    cns_df.update(xx_y_chrom_df)
    cns_df.update(xy_x_chrom_df)
    cns_df.update(xy_y_chrom_df)

    return cns_df


def cns_impute(cns_df, samples_df, method='extend', cn_columns=None, print_info=True):    
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if method ==  'extend':
        return _impute_extend(cns_df, cn_columns, print_info)
    if method == 'diploid':
        return _impute_diploid(cns_df, samples_df, cn_columns, print_info)
    if method == 'zero':
        return fill_nans_with_zeros(cns_df, cn_columns, print_info)
    else:
        msg = f"Unknown imputation method: {method}"
        raise Exception(msg)


def fill_nans_with_zeros(cns_df, cn_columns=None, print_info=True):    
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)    
    res_df = cns_df.copy()
    log_info(print_info, f"Filling {res_df[cn_columns].isna().any(axis=1).sum()} NaN rows with zero")
    # Fully missing chromosomes filled with 0
    for col in cn_columns:
        res_df[col] = res_df[col].fillna(0).astype(int)
    return res_df