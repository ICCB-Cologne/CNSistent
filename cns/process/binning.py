import numpy as np
import pandas as pd
from numba import njit

from cns.process.breakpoints import make_breaks
from cns.utils.conversions import breaks_to_segments, segs_to_chrom_dict
from cns.utils.files import find_cn_cols_if_none
from cns.utils import hg19


def add_cns_loc(cns_df, assembly=hg19):
    cns_df["length"] = (cns_df["end"] - cns_df["start"]).astype(np.uint32)
    cns_df["mid"] = cns_df["start"] + cns_df["length"] // 2
    cns_df["cum_mid"] = cns_df["mid"] + cns_df.apply(lambda x: assembly.chr_starts[x["chrom"]], axis=1)
    return cns_df


def sum_cns(cns_df, cn_columns=None):
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
    return cns_df


def group_bins(cns_df, cn_columns=None, fun_type="mean", assembly=hg19):    
    if fun_type not in ["mean", "max", "min"]:
        raise ValueError("to group bins, fun_type must be one of ['mean', 'max', 'min']")
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if "cum_mid" not in cns_df:
        cns_df = add_cns_loc(cns_df, assembly)
    grouped = cns_df.drop("sample_id", axis=1).groupby(["cum_mid"])
    # calculate mean on grouped except for chrom, where take the first value
    agg_scheme = {
        "chrom": "first",
        "start": "first",
        "end": "first",
        "mid": "first",
        "length": "first" 
    }
    for column in cn_columns:
        agg_scheme[column] = fun_type
    grouped = grouped.agg(agg_scheme).reset_index()
    if len(cn_columns) > 1:
        grouped["total_cn"] = grouped[cn_columns].sum(axis=1)
    return grouped


@njit
def mean_func(cns_array):
    result = []
    for i in range(cns_array.shape[1] - 1):
        result.append(np.average(cns_array[:, i], weights=cns_array[:, -1]))
    return result


@njit
def max_func(cns_array):
    result = []
    for i in range(cns_array.shape[1] - 1):
        result.append(np.max(cns_array[:, i]))
    return result


@njit
def min_func(cns_array):
    result = []
    for i in range(cns_array.shape[1] - 1):
        result.append(np.min(cns_array[:, i]))
    return result


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
def bin_by_segments(cns_df, segments, fun_type="mean", print_progress=True):
    agg_func = _get_agg_func(fun_type)
    chrom_segments = segs_to_chrom_dict(segments)
    cns_df_view = cns_df.set_index(["sample_id", "chrom"])
    new_rows = []
    indices = cns_df_view.index.unique()
    i = 0
    for (sample, chrom), group in cns_df_view.groupby(level=[0, 1]):
        if print_progress:
            i += 1
            print(f"Binning chr ({i}/{len(indices)})", end="\r")
        if chrom in chrom_segments:
            for segment in chrom_segments[chrom]:
                if agg_func != None:
                    bin = _regs_to_bin(sample, chrom, group.values, segment, agg_func)
                    new_rows.append(bin)
                else:
                    bin = _cns_in_seg(sample, chrom, group.values, segment)
                    new_rows.extend(bin)
    if print_progress:
        print("")
    bin_df = pd.DataFrame(new_rows, columns=cns_df.columns)
    bin_df["start"] = bin_df["start"].astype(np.uint32)
    bin_df["end"] = bin_df["end"].astype(np.uint32)
    return bin_df


def bin_by_breaks(cns_df, breaks, fun_type="mean", print_progress=True):
    segments = breaks_to_segments(breaks)
    return bin_by_segments(cns_df, segments, fun_type, print_progress)


def bin_by_break_type(cns_df, break_type, assembly=hg19, fun_type="mean", print_progress=True):
    breaks = make_breaks(break_type, assembly)
    bin_df = bin_by_breaks(cns_df, breaks, fun_type, print_progress)
    return bin_df
