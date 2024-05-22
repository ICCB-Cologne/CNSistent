import numpy as np
import pandas as pd
from numba import njit

from cns.process.breakpoints import get_breakpoints
from cns.utils.conversions import breaks_to_segments
from cns.utils.conversions import segs_to_chrom_dict
from cns.utils import hg19


def add_cns_loc(cns_df, assembly=hg19):
    cns_df = cns_df.copy()
    cns_df["length"] = (cns_df["end"] - cns_df["start"]).astype(np.uint32)
    cns_df["mid"] = cns_df["start"] + cns_df["length"] // 2
    cns_df["cum_mid"] = cns_df["mid"] + cns_df.apply(lambda x: assembly.chr_starts[x["chrom"]], axis=1)
    if "major_cn" and "minor_cn" in cns_df:
        cns_df["total_cn"] = cns_df["major_cn"] + cns_df["minor_cn"]
    if "cn_a" and "cn_b" in cns_df:
        cns_df["total_cn"] = cns_df["cn_a"] + cns_df["cn_b"]
    if "hap_a" and "hap_b" in cns_df:
        cns_df["total_cn"] = cns_df["hap_a"] + cns_df["hap_b"]
    # order by cum_mid
    if "sample_id" in cns_df:
        return cns_df.sort_values(by=["sample_id", "cum_mid"])
    else:
        return cns_df.sort_values(by=["cum_mid"])


# TODO: work with single-column
def mean_bins(bins_df, assembly=hg19):
    if "cum_mid" not in bins_df:
        bins_df = add_cns_loc(bins_df, assembly)
    grouped = bins_df.drop("sample_id", axis=1).groupby(["cum_mid"])
    # calculate mean on grouped except for chrom, where take the first value
    grouped = grouped.agg(
        {
            "chrom": "first",
            "start": "first",
            "end": "first",
            "major_cn": "mean",
            "minor_cn": "mean",
            "mid": "first",
            "length": "first",
        }
    ).reset_index()
    grouped["total_cn"] = grouped["major_cn"] + grouped["minor_cn"]
    return grouped


@njit
def mean_func(major_cn, minor_cn, length):
    return np.average(major_cn, weights=length), np.average(minor_cn, weights=length)


@njit
def max_func(major_cn, minor_cn, length):
    return np.max(major_cn), np.max(minor_cn)


@njit
def min_func(major_cn, minor_cn, length):
    return np.min(major_cn), np.min(minor_cn)


# TODO: work with single-column
def _regs_to_bin(sample_id, chrom, sample_rows, segment, agg_func):
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
        length = end - start
        seg_cns.append((row[2], row[3], length))
        if row[1] >= seg_end:  # last row ends behind the segment
            break
        row_id += 1
        if row_id >= len(sample_rows):
            break
    
    seg_cns = [x for x in seg_cns if not np.isnan(x[0]) and not np.isnan(x[1])]
    if seg_cns == []:
        return (sample_id, chrom, seg_start, seg_end, np.nan, np.nan)
    sel_array = np.array(seg_cns, dtype=np.uint32)
    major_cn, minor_cn = agg_func(sel_array[:, 0], sel_array[:, 1], sel_array[:, 2])
    return (sample_id, chrom, seg_start, seg_end, major_cn, minor_cn)


# TODO: work with single-column
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
        seg_cns.append((sample_id, chrom, start, end, row[2], row[3]))
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
    df_cols = ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"]
    bin_df = pd.DataFrame(new_rows, columns=df_cols)
    bin_df["start"] = bin_df["start"].astype(np.uint32)
    bin_df["end"] = bin_df["end"].astype(np.uint32)
    return bin_df


def bin_by_breaks(cns_df, breaks, fun_type="mean", print_progress=True):
    segments = breaks_to_segments(breaks)
    return bin_by_segments(cns_df, segments, fun_type, print_progress)


def bin_by_break_type(cns_df, break_type, assembly=hg19, fun_type="mean", print_progress=True):
    breaks = get_breakpoints(break_type, assembly)
    bin_df = bin_by_breaks(cns_df, breaks, fun_type, print_progress)
    return bin_df
