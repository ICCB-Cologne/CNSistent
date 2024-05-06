import numpy as np
import pandas as pd
from numba import njit

from cns.process.breakpoints import get_breakpoints
from cns.process.segments import breaks_to_segments
from cns.utils.conversions import segs_to_chrom_dict
from cns.utils.assemblies import hg19


def add_derived(cns_df, assembly=hg19):
    cns_df = cns_df.copy()
    cns_df["total_cn"] = cns_df["major_cn"] + cns_df["minor_cn"]
    cns_df["length"] = (cns_df["end"] - cns_df["start"]).astype(np.uint32)
    cns_df["mid"] = cns_df["start"] + cns_df["length"] // 2
    cns_df["cum_mid"] = cns_df["mid"] + cns_df.apply(lambda x: assembly.chr_starts[x["chrom"]], axis=1)
    # order by cum_mid
    cns_df = cns_df.sort_values(by=["sample_id", "cum_mid"])
    return cns_df


def mean_bins(bins_df, assembly=hg19):
    if "cum_mid" not in bins_df:
        bins_df = add_derived(bins_df, assembly)
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


def _get_bin_for_seg(sample_id, chrom, sample_rows, segment, agg_func):
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


def _get_agg_func(fun_type):
    if fun_type == "mean":
        return mean_func
    if fun_type == "max":
        return max_func
    if fun_type == "min":
        return min_func
    raise ValueError("fun_type must be 'mean', 'max' or 'min', got " + fun_type)


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
                bin = _get_bin_for_seg(sample, chrom, group.values, segment, agg_func)
                new_rows.append(bin)
    if print_progress:
        print("")
    df_cols = ["sample_id", "chrom", "start", "end", "major_cn", "minor_cn"]
    bin_df = pd.DataFrame(new_rows, columns=df_cols)
    return bin_df


def bin_by_breaks(cns_df, breaks, fun_type="mean", print_progress=True):
    segments = breaks_to_segments(breaks)
    return bin_by_segments(cns_df, segments, fun_type, print_progress)


def bin_block(cns_df, break_type, assembly=hg19, fun_type="mean", print_progress=True):
    breaks = get_breakpoints(break_type, assembly)
    bin_df = bin_by_breaks(cns_df, breaks, fun_type, print_progress)
    labelled_bins = add_derived(bin_df, assembly)
    return labelled_bins
