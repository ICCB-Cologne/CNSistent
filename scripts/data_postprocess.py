#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from cns.data_utils import load_samples_out, load_cns_out, out_path
from cns.process.breakpoints import get_breaks  
from cns.process.binning import add_cns_loc, bin_by_breaks, group_bins
from cns.utils.files import get_cn_columns, save_cns, save_samples

def merge_samples(print_info=False):
    # Load the sample data
    prim_samples = load_samples_out("TRACERx_prim_samples.tsv")
    met_samples = load_samples_out("TRACERx_met_samples.tsv")

    # Merge the DataFrames on the common key 'sample_id'
    merged_df = pd.merge(prim_samples, met_samples, on='sample_id', suffixes=('_prim', '_met'))

    # Identify numerical and non-numerical columns
    numerical_cols = merged_df.select_dtypes(include='number').columns
    non_numerical_cols = merged_df.select_dtypes(exclude='number').columns

    # Average the numerical columns
    for col in numerical_cols:
        if col.endswith('_prim'):
            base_col = col[:-5]
            merged_df[base_col] = merged_df[[f'{base_col}_prim', f'{base_col}_met']].mean(axis=1)

    # Take non-numerical values from prim_samples
    for col in non_numerical_cols:
        if col.endswith('_prim'):
            base_col = col[:-5]
            merged_df[base_col] = merged_df[f'{base_col}_prim']

    # Drop the intermediate columns
    merged_df = merged_df[[col for col in merged_df.columns if not col.endswith(('_prim', '_met'))]]

    # Reorder columns to match the order in prim_samples
    merged_df = merged_df[prim_samples.columns]

    # Round down and convert to integer for columns that were integers in prim_samples
    for col in prim_samples.select_dtypes(include='int').columns:
        merged_df[col] = np.floor(merged_df[col]).astype(np.int64)

    # Save the merged DataFrame
    save_samples(merged_df, f"{out_path}/TRACERx_samples.tsv")

def merge_cns(print_info=False):
    prim_cns = load_cns_out("TRACERx_prim_cns_imp.tsv")
    prim_cns.set_index(["sample_id"], inplace=True)
    met_cns = load_cns_out("TRACERx_met_cns_imp.tsv")
    met_cns.set_index(["sample_id"], inplace=True)

    common_samples = list(set(prim_cns.index).intersection(set(met_cns.index)))
    filtered_prim_cns = prim_cns[~prim_cns.index.isin(common_samples)].reset_index()
    filtered_met_cns = met_cns[~met_cns.index.isin(common_samples)].reset_index()
    binned_samples = [filtered_prim_cns, filtered_met_cns]

    i = 1
    for sample_id in common_samples:
        if print_info:
            print(f"Merging for sample {i}/{len(common_samples)}", end="\r")
        i = i + 1
        sel_prim = prim_cns.loc[sample_id].copy().reset_index()
        sel_met = met_cns.loc[sample_id].copy().reset_index()
        prim_breaks = get_breaks(sel_prim)
        met_breaks = get_breaks(sel_met)
        merged_breaks = { chrom: np.sort(np.unique(prim_breaks[chrom] + met_breaks[chrom])) for chrom in prim_breaks.keys() }
        binned_prim = add_cns_loc(bin_by_breaks(sel_prim, merged_breaks, print_progress=False))
        binned_met = add_cns_loc(bin_by_breaks(sel_met, merged_breaks, print_progress=False))
        binned_merge = group_bins(pd.concat([binned_prim, binned_met]))
        cn_columns = get_cn_columns(binned_merge)
        binned_merge["sample_id"] = sample_id
        binned_merge = binned_merge[["sample_id", "chrom", "start", "end"] + cn_columns]
        binned_samples.append(binned_merge)

    binned_df = pd.concat(binned_samples).sort_values(["sample_id", "chrom", "start"])
    save_cns(binned_df, f"{out_path}/TRACERx_cns_imp.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts source datasets into format usable by cns")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print debug information")
    args = parser.parse_args()
    print_debug = args.verbose
    if print_debug:
        print("Merging samples for TRACERx...")
    merge_samples(print_debug)
    if print_debug:
        print("Merging CNS for TRACERx...")
    merge_cns(print_debug)