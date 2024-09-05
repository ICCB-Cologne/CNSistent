import numpy as np
import pandas as pd
from numba import njit

from cns.process.breakpoints import make_breaks
from cns.utils.conversions import breaks_to_segments, segs_to_chrom_dict
from cns.utils.files import find_cn_cols_if_none, get_cn_columns
from cns.utils import hg19
from cns.utils.logging import log_info


def add_cns_loc(cns_df, assembly=hg19):
    cns_df["length"] = (cns_df["end"] - cns_df["start"]).astype(np.uint32)
    cns_df["mid"] = cns_df["start"] + cns_df["length"] // 2
    offset = cns_df.apply(lambda x: assembly.chr_starts[x["chrom"]], axis=1)
    cns_df["cum_mid"] = cns_df["mid"] + offset
    return cns_df


def sum_cns(cns_df, cn_columns=None):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    # remove total_cn from cn_columns if it is there
    if "total_cn" in cn_columns:
        cn_columns.remove("total_cn")
    cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
    return cns_df


# TODO: Add empty check
def group_bins(cns_df, cn_columns=None, fun_type="mean", group_col='cum_mid', assembly=hg19):    
    if fun_type not in ["mean", "max", "min"]:
        raise ValueError("to group bins, fun_type must be one of ['mean', 'max', 'min']")
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if group_col not in cns_df.columns:
        raise ValueError(f"cns_df must have a column '{group_col}' to group the bins on")
    grouped = cns_df.drop("sample_id", axis=1).groupby([group_col])

    # calculate mean on grouped except for chrom, where take the first value
    agg_scheme = {
        "chrom": "first",
        "start": "first",
        "end": "first",
    }
    if "mid" in cns_df.columns:
        agg_scheme["mid"] = "first"
    if "length" in cns_df.columns:
        agg_scheme["length"] = "first"
    for column in cn_columns:
        agg_scheme[column] = fun_type
    grouped = grouped.agg(agg_scheme).reset_index()
    return grouped


@njit
def mean_func(cns_array):
    return [np.average(cns_array[:, i], weights=cns_array[:, -1]) for i in range(cns_array.shape[1] - 1)]

@njit
def max_func(cns_array):
    return [np.max(cns_array[:, i]) for i in range(cns_array.shape[1] - 1)]

@njit
def min_func(cns_array):
    return [np.min(cns_array[:, i]) for i in range(cns_array.shape[1] - 1)]

def _regs_to_bin(sample_id, chrom, sample_rows, segment, agg_func):
    row_id = 0
    seg_cns = []
    cns_cols = sample_rows.shape[1] - 2
    seg_start, seg_end = segment
    while sample_rows[row_id][1] <= seg_start:
        row_id += 1
        if row_id >= len(sample_rows):
            break
    while sample_rows[row_id][0] < seg_end:
        row = sample_rows[row_id]
        start = max(row[0], seg_start)
        end = min(row[1], seg_end)
        seg_cns.append(np.concatenate([row[2:],  [end - start]]))
        if row[1] >= seg_end:  # last row ends behind the segment
            break
        row_id += 1
        if row_id >= len(sample_rows):
            break
    
    seg_cns = [x for x in seg_cns if not np.isnan(x[0]) and not np.isnan(x[1])]
    
    if seg_cns == []:
        return [sample_id, chrom, seg_start, seg_end] + [np.nan] * cns_cols
    sel_array = np.array(seg_cns, dtype=np.uint32)
    cns = agg_func(sel_array)
    return [sample_id, chrom, seg_start, seg_end] + cns


def _cns_in_seg(sample_id, chrom, sample_rows, segment):
    row_id = 0
    seg_cns = []
    seg_start, seg_end = segment
    while sample_rows[row_id][1] <= seg_start:
        row_id += 1
        if row_id >= len(sample_rows):
            break
    while sample_rows[row_id][0] < seg_end:
        row = sample_rows[row_id]
        start = max(row[0], seg_start)
        end = min(row[1], seg_end)
        seg_cns.append([sample_id, chrom, start, end] + list(row[2:]))
        if row[1] >= seg_end:  # last row ends behind the segment
            break
        row_id += 1
        if row_id >= len(sample_rows):
            break

    return seg_cns


def _get_agg_func(fun_type):
    if fun_type == "" or fun_type is None or fun_type == "none":
        return None
    if fun_type == "mean":
        return mean_func
    if fun_type == "max":
        return max_func
    if fun_type == "min":
        return min_func
    raise ValueError("fun_type must be one of ['mean', 'max', 'min', 'none', '']  got " + fun_type)


# Add column names
def bin_by_segments(cns_df, segments, fun_type="mean", print_info=True):
    agg_func = _get_agg_func(fun_type)
    chrom_segments = segs_to_chrom_dict(segments)
    cn_columns = get_cn_columns(cns_df)
    sel_cols = ["sample_id", "chrom", "start", "end"] + cn_columns
    cns_df_view = cns_df[sel_cols].set_index(["sample_id", "chrom"])
    new_rows = []
    indices = cns_df_view.index.unique()
    i = 0
    for i, ((sample, chrom), group) in enumerate(cns_df_view.groupby(level=[0, 1])):
        log_info(print_info, print(f"Binning chr ({i+1}/{len(indices)})", end="\r"))
        if chrom in chrom_segments:
            for segment in chrom_segments[chrom]:
                if agg_func != None:
                    bin = _regs_to_bin(sample, chrom, group.values, segment, agg_func)
                    new_rows.append(bin)
                else:
                    bin = _cns_in_seg(sample, chrom, group.values, segment)
                    new_rows.extend(bin)
    log_info(print_info, print(f"Binning finished. Converting {len(new_rows)} rows...", end="\r"))
    bin_df = pd.DataFrame(new_rows, columns=sel_cols)
    bin_df["start"] = bin_df["start"].astype(np.uint32)
    bin_df["end"] = bin_df["end"].astype(np.uint32)
    log_info(print_info, print(f"Binned into {len(new_rows)} CNS." + " " * 40))
    return bin_df


def bin_by_breaks(cns_df, breaks, fun_type="mean", print_info=True):
    segments = breaks_to_segments(breaks)
    return bin_by_segments(cns_df, segments, fun_type, print_info)


def bin_by_break_type(cns_df, break_type, assembly=hg19, fun_type="mean", print_info=True):
    breaks = make_breaks(break_type, assembly)
    bin_df = bin_by_breaks(cns_df, breaks, fun_type, print_info)
    return bin_df
