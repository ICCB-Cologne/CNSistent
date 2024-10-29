#!/usr/bin/env python

from multiprocessing import Pool
import pandas as pd
import time
import argparse
from os.path import exists

from cns.process.segments import regions_remove, regions_select
from cns.utils.assemblies import get_assembly
from cns.utils.canonization import find_cn_cols_if_none
from cns.utils.files import load_cns, load_segments, save_cns, save_segments, dataframe_array_split, samples_df_from_cns_df, load_samples, fill_sex_if_missing, save_samples
from cns.process.pipelines import main_fill, main_impute, main_aggregate, main_coverage, main_ploidy, main_segment, main_signatures
from cns.utils.logging import log_info, log_warn


def _add_sp_args(action, parser):
    parser.add_argument("data", type=str, help="path to the TSV file with copy number segments")    
    parser.add_argument("--samples", type=str, help="path to the samples file", required=False, default="")
    parser.add_argument("--out", type=str, help="output file path", required=False, default="./cns.out.tsv")
    parser.add_argument("--assembly", type=str, help="assembly to use. One of: hg19, hg38.", required=False, default="hg19")
    parser.add_argument("--threads", type=int, help="number of threads to use", required=False, default=1)
    parser.add_argument("--verbose", help="print progress to console", action="store_true")
    parser.add_argument("--time", help="save runtime info", action="store_true")
    parser.add_argument("--subsplit", type=int, help="will split the processing into chunks to lower memory needs", required=False, default=1)
    parser.add_argument("--cncols", type=str, help="The name of either a single CN column or two comma separated columns.", required=False, default=None)

    if action in ["coverage", "ploidy", "signatures"]:
        parser.add_argument("--segments", type=str, help="A path to a segmentation file defining the segments to compute on.", required=False, default=None)

    if action == "aggregate":
        parser.add_argument("--segments", type=str, help="A path to a segmentation file defining the segments to compute on.", required=True)
        parser.add_argument("--how", type=str, help="The aggregation function, one of ['min', 'max', 'mean']", required=False, default="mean")

    if action == "segment":
        parser.add_argument("--split", type=int, help="Distance in which regions should be, can be a positive integer or 0 for no splitting (whole regions).", required=False, default=0)
        parser.add_argument("--select", type=str, help="Selects the regions to create segments based  on, can be either 'arms', 'bands', path to a BED file, or empty for whole chromosomes.", required=False, default="")
        parser.add_argument("--remove", type=str, help="Removed the regions after selection, before segmentation, can be either 'gaps', path to a BED file, or empty.", required=False, default="")
        parser.add_argument("--filter", type=int, help="If set, regions smaller than the given size are excluded from selection and gaps.", required=False, default=0)
        parser.add_argument("--merge", type=int, help="Maximum distance between breakpoint clusters for breakpoint merging", required=False, default=0)

def _parse_args():
    # Top-level parser
    parser = argparse.ArgumentParser(description="Impute missing values in CNS data")
    subparsers = parser.add_subparsers(dest="action", help="cns action to perform. One of:")

    parser.add_argument("-v", action="version", version="%(prog)s 1.0")

    sp_dict = {}
    sp_dict["fill"] = subparsers.add_parser("fill", help=f"Adds Nan regions to the CNS data to match the assembly.")
    sp_dict["impute"] = subparsers.add_parser("impute", help=f"Imputes missing values in the CNS data.")
    sp_dict["coverage"] = subparsers.add_parser("coverage", help=f"Calculates coverage for filled (but not imputed) CNS data.")
    sp_dict["ploidy"] = subparsers.add_parser("ploidy", help=f"Conducts breakpoint analysis for CNS data (NaNs are ignored).")
    sp_dict["signatures"] = subparsers.add_parser("signatures", help=f"Extracts basal CN signatures from CNS data (NaNs are ignored).")
    sp_dict["segment"] = subparsers.add_parser("segment", help=f"Calculates segmentation regions for CNS data.")
    sp_dict["aggregate"] = subparsers.add_parser("aggregate", help=f"Aggregate copy numbers across segments to fill provided segments.")    
    for action, sp in sp_dict.items():
        _add_sp_args(action=action, parser=sp)

    args = parser.parse_args()
    if args.action is None:
        parser.print_help()
        exit(1)
    if args.action not in sp_dict:
        raise ValueError(f"Action {args.action} not recognized.")
    
    if args.action == "segment":
        if args.merge > 0 and not args.data:
            raise ValueError("Merging breakpoints requires the --data option to provide CNS data.")
        if args.split < 0:
            args.split = 0
        if args.select not in ["", "arms", "bands"] and not exists(args.select):
            raise ValueError(f"Selection {args.select} is not a build-in or a path to a file.")
        if args.remove not in ["", "gaps"] and not exists(args.remove):
            raise ValueError(f"Removal {args.remove} is not a build-in or a path to a file.")
        if args.threads > 1:
            print("segmentation is not data parallelizable, --threads option will be ignored.")
            args.threads = 1
        if args.subsplit > 1:
            print("segmentation is not data parallelizable, --subsplit option will be ignored.")
            args.subsplit = 1            
    
    if args.threads <= 0:
        raise ValueError("The --threads option must be greater than 0.")
    
    if args.subsplit <= 0:
        raise ValueError("The --subsplit option must be greater than 0.")

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
    elif action == "signatures":
        return main_signatures  
    elif action == "segment":
        return main_segment
    elif action == "aggregate":
        return main_aggregate
    else:
        raise ValueError(f"Action {action} not recognized.")


def _get_segs_df(segs_arg):
    if segs_arg != "" and segs_arg != None:
        return load_segments(segs_arg)
    return None


def _get_blocks(action, cns_blocks, samples_blocks, cols_block, assembly, args):
    block_count = len(cns_blocks)
    # Apply process_block to each pair of blocks
    ass_block = [assembly]*block_count
    ver_block = [False]*block_count
    ver_block[-1] = args.verbose
    cols_block = [cols_block]*block_count
    if action == "impute":
        ext_block = ['extend']*block_count
        return zip(cns_blocks, samples_blocks, ext_block, cols_block, ver_block)
    if action == "fill":
        add_missing = [True]*block_count
        return zip(cns_blocks, samples_blocks, cols_block, ass_block, add_missing, ver_block)        
    elif action == "segment":
        select = regions_select(args.select, assembly)
        remove = regions_remove(args.remove, assembly)
        select_block = [select]*block_count
        remove_block = [remove]*block_count
        split_block = [args.split]*block_count
        merge_block = [args.merge]*block_count
        filter_block = [args.filter]*block_count
        return zip(cns_blocks, select_block, remove_block, split_block, merge_block, filter_block, ass_block, ver_block)
    elif action in ["coverage", "ploidy", "signatures"]:
        segs_block = [_get_segs_df(args.segments)] * block_count
        return zip(cns_blocks, samples_blocks, cols_block, segs_block, ass_block, ver_block)
    elif action == "aggregate":
        segs_block = [_get_segs_df(args.segments)] * block_count
        fun_block = [args.how]*block_count
        return zip(cns_blocks, segs_block, fun_block, cols_block, ver_block)
    else:
        raise ValueError(f"Unknown action {action}")


def _process(action, cns_df, samples_df, cn_cols, assembly, args):   
    main_fun = _action_to_fun(action)
    threads = abs(args.threads)
    samples_blocks = dataframe_array_split(samples_df, threads)
    cns_blocks = [cns_df.query("sample_id in @block.index").reset_index(drop=True) for block in samples_blocks]
    zip_blocks = _get_blocks(action, cns_blocks, samples_blocks, cn_cols, assembly, args)
    if threads == 1:
        return [main_fun(*list(*zip_blocks))]
    else:
        with Pool(threads) as pool:
            log_info(args.verbose, f"Multiprocessing with {threads} threads..")
            res_blocs = pool.starmap(main_fun, zip_blocks)        
            pool.close()   
            pool.join()    
        return res_blocs
    

def parse_cncols(cncols):
    if cncols != None:
        cncols = cncols.split(",")
        if len(cncols) > 2:
            raise ValueError("Only one or two columns can be specified.")
    return cncols


def main():
    args = _parse_args()
    action = args.action
    assembly = get_assembly(args.assembly)
    cns_file_path = args.data
    samples_path = args.samples
    out_file = args.out
    print_info = args.verbose
    subsplit = args.subsplit
    cncols = parse_cncols(args.cncols)

    # Read the input
    log_info(print_info, f"***** cns {action} *****")

    # For segmentation without cns file we don't use cns files and multiprocessing
    if action == "segment" and args.merge == 0:
        select = regions_select(args.select, assembly)
        remove = regions_remove(args.remove, assembly)
        segs = main_segment(None, select, remove, args.split, args.merge, args.filter, assembly, print_info)        
        save_segments(segs, out_file)
        log_info(print_info, "Done.")
        return

    if not exists(cns_file_path):    
        raise ValueError(f"Copy number file {cns_file_path} not found.")
    else:
        log_info(print_info, f"CNS file at {cns_file_path}...")     
    cns_df = load_cns(cns_file_path, canonize=True, cn_columns=cncols, assembly=assembly, print_info=print_info)
    cn_columns = find_cn_cols_if_none(cns_df, cncols)     
    if samples_path == "":
        samples_df = samples_df_from_cns_df(cns_df)
    else:
        samples_df = load_samples(samples_path)
        samples_df = fill_sex_if_missing(cns_df, samples_df)
    samples_blocks = dataframe_array_split(samples_df, subsplit) 

    for i in range(subsplit):
        log_info(print_info, f"Processing block {i+1}/{subsplit}...")
        
        samples_block = samples_blocks[i]
        
        # Perform the action
        start = time.time()
        res_list = _process(action, cns_df, samples_block, cn_columns, assembly, args)
        runtime = time.time() - start

        if print_info:
            print(f"Finished in {runtime:.3f} seconds. Writing to {out_file}...")        
            if args.time:  
                timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))  
                with open("./out/times.tsv", "a") as f:
                    f.write(f"{timestamp}\t{action}\t{args.threads}\t{out_file}\t{runtime}\n")

        for j in range(len(res_list)): 
            mode = "w" if i == 0 and j == 0 else "a"
            res = res_list[j]
            if action == "segment":
                print(out_file)
                save_segments(res, out_file)
            elif action in ["fill", "impute", "aggregate"]:
                save_cns(res, out_file, change_coords=True, mode=mode)
            elif action in ["coverage", "ploidy", "signatures"]:
                save_samples(res, out_file, mode=mode)
            else:
                raise ValueError(f"Unknown action {action}")
    
    log_info(print_info, "Done.")


if __name__ == "__main__":
    main()
