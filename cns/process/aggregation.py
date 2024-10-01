import numpy as np
import pandas as pd
from numba import njit

from cns.process.breakpoints import make_breaks
from cns.utils.conversions import breaks_to_segments
from cns.utils.canonization import find_cn_cols_if_none
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
# TODO: Group on start, end
def group_samples(cns_df, cn_columns=None, how="mean", group_col='cum_mid', assembly=hg19):    
    if how not in ["mean", "max", "min"]:
        raise ValueError("to group samples, how must be one of ['mean', 'max', 'min']")
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    if group_col not in cns_df.columns:
        raise ValueError(f"cns_df must have a column '{group_col}' to group the samples on")
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
        agg_scheme[column] = how
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

def _aggregate_regs(sample_id, chrom, sample_rows, seg_start, seg_end, seg_name, agg_func):
    row_id = 0
    seg_cns = []
    cns_cols = sample_rows.shape[1] - 2
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
    return [sample_id, chrom, seg_start, seg_end] + cns + [seg_name]


def _mask_by_regs(sample_id, chrom, sample_rows, seg_start, seg_end, seg_name):
    row_id = 0
    seg_cns = []
    while sample_rows[row_id][1] <= seg_start:
        row_id += 1
        if row_id >= len(sample_rows):
            break
    while row_id < len(sample_rows) and sample_rows[row_id][0] < seg_end:
        row = sample_rows[row_id]
        start = max(row[0], seg_start)
        end = min(row[1], seg_end)
        seg_cns.append([sample_id, chrom, start, end] + list(row[2:]) + [seg_name])
        if row[1] >= seg_end:  # last row ends behind the segment
            break
        row_id += 1
        if row_id >= len(sample_rows):
            break
    return seg_cns


def _get_agg_func(how):
    if how == "" or how is None or how == "none":
        return None
    if how == "mean":
        return mean_func
    if how == "max":
        return max_func
    if how == "min":
        return min_func
    raise ValueError("how must be one of ['mean', 'max', 'min', 'none', '']  got " + how)


# Add column names
def aggregate_by_segments(cns_df, segs, how="mean", cn_columns=None, print_info=True):
    agg_func = _get_agg_func(how)
    cn_columns = find_cn_cols_if_none(cns_df, cn_columns)
    sel_cols = ["sample_id", "chrom", "start", "end"] + cn_columns
    cns_df_view = cns_df[sel_cols].set_index(["sample_id", "chrom"])
    new_rows = []
    indices = cns_df_view.index.unique()
    i = 0
    for i, ((sample, chrom), group) in enumerate(cns_df_view.groupby(level=[0, 1])):
        if print_info:
            print(f"Aggregating chr ({i+1}/{len(indices)})", end="\r")
        if chrom not in segs:
            continue
        for seg_start, seg_end, seg_name in segs[chrom]:
            if agg_func != None:
                cn_segs = _aggregate_regs(sample, chrom, group.values, seg_start, seg_end, seg_name, agg_func)
                new_rows.append(cn_segs)
            else:
                cn_segs = _mask_by_regs(sample, chrom, group.values, seg_start, seg_end, seg_name)
                new_rows.extend(cn_segs)
    if print_info:
        print(f"Aggregation finished. Converting {len(new_rows)} rows...", end="\r")
    res_df = pd.DataFrame(new_rows, columns=sel_cols + ["name"]) 
    res_df["start"] = res_df["start"].astype(np.uint32)
    res_df["end"] = res_df["end"].astype(np.uint32)
    log_info(print_info, f"Aggregated into {len(new_rows)} CNS." + " " * 40)
    return res_df


def aggregate_by_breaks(cns_df, breaks, how="mean", cn_columns=None, print_info=True):
    segments = breaks_to_segments(breaks)
    return aggregate_by_segments(cns_df, segments, how, cn_columns, print_info)


def aggregate_by_break_type(cns_df, break_type, assembly=hg19, how="mean", cn_columns=None, print_info=True):
    breaks = make_breaks(break_type, assembly)
    return aggregate_by_breaks(cns_df, breaks, how, cn_columns, print_info)
