#!/usr/bin/env python

from multiprocessing import Pool
import numpy as np
import pandas as pd
import time
import argparse
from os.path import join as exists

from cns.process.binning import bin_by_segments
from cns.process.imputation import fill_sex_if_missing
from cns.process.segments import get_genome_segments
from cns.utils.assemblies import get_assembly
from cns.utils.files import load_cns, save_cns, save_regions
from cns.process.pipelines import main_fill, main_impute, main_coverage, main_ploidy, main_cluster, regions_remove, regions_select
from cns.utils.files import dataframe_array_split, samples_df_from_cns_df, load_samples


def _add_common_args(parser):
    parser.add_argument("data", type=str, help="Path to the TSV file with copy number segments")
    parser.add_argument("--samples", type=str, help="Path to the samples file", required=False, default="")
    parser.add_argument("--out", type=str, help="Output file.", required=False, default="./cns.out.tsv")
    parser.add_argument("--assembly", type=str, help="Assembly to use. One of: hg19, hg38.", required=False, default="hg19")
    parser.add_argument("--threads", type=int, help="Number of threads to use", required=False, default=1)
    parser.add_argument("--verbose", help="Print progress", action="store_true")
    parser.add_argument("--time", help="Save runtime", action="store_true")


def _parse_args():
    # Top-level parser
    parser = argparse.ArgumentParser(description="Impute missing values in CNS data")
    subparsers = parser.add_subparsers(dest="action", help="cns action to perform. One of:")

    parser.add_argument("-v", action="version", version="%(prog)s 1.0")

    sp_dict = {}
    sp_dict["fill"] = subparsers.add_parser("fill", help=f"Adds Nan regions to the CNS data to match the assembly.")
    sp_dict["impute"] = subparsers.add_parser("impute", help=f"Imputes missing values in the CNS data.")
    sp_dict["coverage"] = subparsers.add_parser("coverage", help=f"Calculates coverage for filled (but not imputed) CNS data.")
    sp_dict["ploidy"] = subparsers.add_parser("ploidy", help=f"Calculates aneuploidy for CNS data (NaNs are ignored).")
    sp_dict["bin"] = subparsers.add_parser("bin", help=f"Creates bins for CNS data.")
    sp_dict["cluster"] = subparsers.add_parser("cluster", help=f"Calculates a clustering of breakpoints between samples.")

    for sp in sp_dict.values():
        _add_common_args(sp)

    # TODO: Add argument to remove a file other than gaps
    # TODO: Move arms/bands to regions
    sp_dict["bin"].add_argument("--bins", type=int, help="Size of the bins, can be a positive integer or 0 for no bins (whole regions).", required=False, default=0)
    sp_dict["bin"].add_argument("--select", type=str, help="Selects the regions to bin on, can be either 'arms', 'bands', path to a BED file, or empty for whole chromosomes.", required=False, default="")
    sp_dict["bin"].add_argument("--remove", type=str, help="Removed the regions after selection, before binning, can be either 'gaps', path to a BED file, or empty.", required=False, default="")
    sp_dict["bin"].add_argument("--filter", type=int, help="If set, regions smaller than the given size are exclued from selection and gaps.", required=False, default=0)
    sp_dict["bin"].add_argument("--aggregate", type=str, help="The aggregation function, one of ['min', 'max', 'mean']", required=False, default="mean")
    sp_dict["bin"].add_argument("--onlybins", help="If set, creates a BED file with regions corresponding to the individual bins, without the actual binning. (A dummy CNS file still has to be provided)", action="store_true")
    
    sp_dict["cluster"].add_argument("--dist", type=int, help="Maximum distance between breakpoint clusters for merging", required=False, default=0)

    args = parser.parse_args()
    if args.action is None:
        parser.print_help()
        exit(1)
    if args.action not in sp_dict:
        raise ValueError(f"Action {args.action} not recognized.")
    
    if args.action == "bin":
        if args.bins < 0:
            args.bins = 0
        if args.select not in ["", "arms", "bands"] and not exists(args.select):
            raise ValueError(f"Selection {args.select} is not a build-in or a path to a file.")
        if args.remove not in ["", "gaps"] and not exists(args.remove):
            raise ValueError(f"Removal {args.remove} is not a build-in or a path to a file.")

    if args.action == "cluster":
        if args.threads > 1:
            print("WARNING: Clustering is not data parallelizable, --threads option will be ignored.")
            args.threads = 1

    return args


def _action_to_fun(action):
    if action == "fill":
        return main_fill
    elif action == "impute":
        return main_impute
    elif action == "coverage":
        return main_coverage
    elif action == "ploidy":
        return main_ploidy   
    elif action == "bin":
        return bin_by_segments
    elif action == "cluster":
        return main_cluster
    else:
        raise ValueError(f"Action {action} not recognized.")


def _get_segments(args, assembly):    
    bin_size = args.bins
    filter_size = args.filter
    select = regions_select(args.select, assembly)
    remove = regions_remove(args.remove, assembly)
    segs = get_genome_segments(select, bin_size, remove, filter_size)
    return segs


def _get_blocks(action, cns_blocks, samples_blocks, assembly, threads, args):
    # Apply process_block to each pair of blocks
    ass_block = [assembly]*threads
    ver_block = [False]*threads
    ver_block[0] = args.verbose
    if action in ["impute"]:
        column_block = [['major_cn', 'minor_cn']]*threads
        return zip(cns_blocks, column_block, ver_block)
    if action in ["fill"]:
        column_block = [['major_cn', 'minor_cn']]*threads
        add_missing = [True]*threads
        return zip(cns_blocks, samples_blocks, ass_block, column_block, add_missing, ver_block)        
    elif action in ["cluster"]:
        dist_block = [args.dist]*threads
        return zip(cns_blocks, dist_block, ass_block, ver_block)
    elif action in ["coverage", "ploidy"]:
        return zip(cns_blocks, samples_blocks, ass_block, ver_block)
    elif action == "bin":
        segs = _get_segments(args, assembly)
        segs_block = [segs]*threads
        fun_block = [args.aggregate]*threads
        return zip(cns_blocks, segs_block, fun_block, ver_block)
    else:
        raise ValueError(f"Unknown action {action}")


def _process(action, cns_df, samples_df, assembly, args):   
    main_fun = _action_to_fun(action)
    threads = args.threads
    samples_blocks = dataframe_array_split(samples_df, threads)
    cns_blocks = [cns_df.query("sample_id in @block.index").reset_index(drop=True) for block in samples_blocks]
    if threads == 1:
        zip_blocks = _get_blocks(action, cns_blocks, samples_blocks, assembly, threads, args)
        return main_fun(*list(*zip_blocks))
    else:
        with Pool(threads) as pool:
            zip_blocks = _get_blocks(action, cns_blocks, samples_blocks, assembly, threads, args)
            res_blocs = pool.starmap(main_fun, zip_blocks)            
        return pd.concat(res_blocs)
    

def main():
    args = _parse_args()
    action = args.action
    assembly = get_assembly(args.assembly)
    cns_file_path = args.data
    samples_path = args.samples
    out_file = args.out
    print_progress = args.verbose

    # Read the input
    if print_progress:
        print(f"cns {action}.")
        print("Parsing input...")

    if not exists(cns_file_path):
        raise ValueError(f"File {cns_file_path} not found.")
    cns_df = load_cns(cns_file_path)

    if samples_path == "":
        samples_df = samples_df_from_cns_df(cns_df)
    else:
        samples_df = load_samples(samples_path)
        samples_df = fill_sex_if_missing(cns_df, samples_df)

    # Perform the action
    if action == "bin" and args.onlybins:        
        if print_progress:
            print(f"Calculating binning segments...")  
        segs = _get_segments(args, assembly)
        res_df = pd.DataFrame(segs, columns=["chrom", "start", "end"])
    else:
        if print_progress:
            print(f"Processing {cns_file_path}...")    
        start = time.time()
        res_df = _process(action, cns_df, samples_df, assembly, args)

    # write out the results
    end = time.time()
    if print_progress:
        runtime = end - start
        print(f"Finished in {runtime:.3f} seconds.")
        print(f"Writing to {out_file}...")
        # Add to file times.tsv
        if args.time:
            with open("./out/times.tsv", "a") as f:
                f.write(f"{action}\t{args.threads}\t{out_file}\t{runtime}\n")
        
    if action == "bin" and args.onlybins:
        save_regions(res_df, out_file, True)
    elif action in ["fill", "impute", "bin"]:
        save_cns(res_df, out_file, sort=True)
    else:
        res_df.sort_index()
        res_df.to_csv(out_file, sep="\t", index=True)

    if print_progress:
        print("Done.")


if __name__ == "__main__":
    main()