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

    # Filtered indices
    common_samples = merged_df.index.unique()
    met_only = met_samples[~met_samples.index.isin(common_samples)].copy()

    # dataset label
    prim_samples["dataset"] = "primary"
    prim_samples[prim_samples.isin(common_samples)]["dataset"] = "both"	
    met_only["dataset"] = "metastatic"

    all_df = pd.concat([prim_samples,  met_only], axis=0)

    # Save the merged DataFrame
    save_samples(all_df, f"{out_path}/TRACERx_samples.tsv")


def merge_cns(print_info=False, filled=False):
    suffix = "fill" if filled else "imp"
    prim_cns = load_cns_out(f"TRACERx_prim_cns_{suffix}.tsv")
    prim_cns.set_index(["sample_id"], inplace=True)
    met_cns = load_cns_out(f"TRACERx_met_cns_{suffix}.tsv")
    met_cns.set_index(["sample_id"], inplace=True)

    common_samples = list(set(prim_cns.index).intersection(set(met_cns.index)))
    filtered_met_cns = met_cns[~met_cns.index.isin(common_samples)].reset_index()
    prim_cns.reset_index(inplace=True)

    binned_df = pd.concat([prim_cns, filtered_met_cns]).sort_values(["sample_id", "chrom", "start"])
    save_cns(binned_df, f"{out_path}/TRACERx_cns_{suffix}.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts source datasets into format usable by cns")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print debug information")
    args = parser.parse_args()
    print_debug = args.verbose
    if print_debug:
        print("Merging samples for TRACERx...")
    merge_samples(print_debug)
    if print_debug:
        print("Merging filled CNS for TRACERx...")
    merge_cns(print_debug, True)    
    if print_debug:
        print("Merging imputed CNS for TRACERx...")
    merge_cns(print_debug, False)