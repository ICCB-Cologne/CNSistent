#!/usr/bin/env python

from multiprocessing import Pool
import pandas as pd
import time
import argparse
from os.path import exists

from cns.process.segments import regions_remove, regions_select
from cns.utils.assemblies import get_assembly
from cns.utils.canonization import find_cn_cols_if_none
from cns.utils.files import load_cns, load_regions, save_cns, save_regions, dataframe_array_split, samples_df_from_cns_df, load_samples, fill_sex_if_missing
from cns.process.pipelines import main_fill, main_impute, main_bin, main_coverage, main_ploidy, main_segment
from cns.utils.logging import log_info


def _add_common_args(parser):
    parser.add_argument("data", type=str, help="path to the TSV file with copy number segments") #TODO: make optional
    parser.add_argument("--samples", type=str, help="path to the samples file", required=False, default="")
    parser.add_argument("--out", type=str, help="output file path", required=False, default="./cns.out.tsv")
    parser.add_argument("--assembly", type=str, help="assembly to use. One of: hg19, hg38.", required=False, default="hg19")
    parser.add_argument("--threads", type=int, help="number of threads to use", required=False, default=1)
    parser.add_argument("--verbose", help="print progress to console", action="store_true")
    parser.add_argument("--time", help="save runtime info", action="store_true")
    parser.add_argument("--noheader", help="if set, the input/output file does not have a header", action="store_true")
    parser.add_argument("--nosample", help="if set, the input/output file does not have an sample column. Not compatible with the --samples option.", action="store_true")
    parser.add_argument("--subsplit", type=int, help="will split the processing into chunks to lower memory needs", required=False, default=1)
    parser.add_argument("--cncols", type=str, help="The name of either a single CN column or two comma separated columns.", required=False, default=None)


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
    sp_dict["segment"] = subparsers.add_parser("cluster", help=f"Calculates segmentation regions for CNS data.")

    for sp in sp_dict.values():
        _add_common_args(sp)

    sp_dict["segment"].add_argument("--split", type=int, help="Distance in which regions should be, can be a positive integer or 0 for no splitting (whole regions).", required=False, default=0)
    sp_dict["segment"].add_argument("--select", type=str, help="Selects the regions to bin on, can be either 'arms', 'bands', path to a BED file, or empty for whole chromosomes.", required=False, default="")
    sp_dict["segment"].add_argument("--remove", type=str, help="Removed the regions after selection, before binning, can be either 'gaps', path to a BED file, or empty.", required=False, default="")
    sp_dict["segment"].add_argument("--filter", type=int, help="If set, regions smaller than the given size are exclued from selection and gaps.", required=False, default=0)
    sp_dict["segment"].add_argument("--merge", type=int, help="Maximum distance between breakpoint clusters for breakpoint merging", required=False, default=0)

    sp_dict["bin"].add_argument("--segments", type=str, help="A path to a segmentation file defining the regions to bin into.", required=True)
    sp_dict["bin"].add_argument("--aggregate", type=str, help="The aggregation function, one of ['min', 'max', 'mean']", required=False, default="mean")

    args = parser.parse_args()
    if args.action is None:
        parser.print_help()
        exit(1)
    if args.action not in sp_dict:
        raise ValueError(f"Action {args.action} not recognized.")
    
    if args.action == "segment":
        if args.split < 0:
            args.split = 0
        if args.select not in ["", "arms", "bands"] and not exists(args.select):
            raise ValueError(f"Selection {args.select} is not a build-in or a path to a file.")
        if args.remove not in ["", "gaps"] and not exists(args.remove):
            raise ValueError(f"Removal {args.remove} is not a build-in or a path to a file.")

    if args.action == "segment":
        if args.threads > 1:
            print("segmentation is not data parallelizable, --threads option will be ignored.")
            args.threads = 1
        if args.subsplit > 1:
            print("segmentation is not data parallelizable, --subsplit option will be ignored.")
            args.subsplit = 1        

    if args.nosample and args["samples"] != "":
        raise ValueError("The --nosample and --samples options are incompatible.")
    
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
    elif action == "bin":
        return main_bin
    elif action == "segment":
        return main_segment
    else:
        raise ValueError(f"Action {action} not recognized.")


def _get_blocks(action, cns_blocks, samples_blocks, cols_block, assembly, threads, args):
    # Apply process_block to each pair of blocks
    ass_block = [assembly]*threads
    ver_block = [False]*threads
    ver_block[-1] = args.verbose
    cols_block = [cols_block]*threads
    if action == "impute":
        ext_block = ['extend']*threads
        return zip(cns_blocks, samples_blocks, ext_block, cols_block, ver_block)
    if action == "fill":
        add_missing = [True]*threads
        return zip(cns_blocks, samples_blocks, cols_block, ass_block, add_missing, ver_block)        
    elif action == "segment":
        select = regions_select(args.select, assembly)
        remove = regions_remove(args.remove, assembly)
        select_block = [select]*threads
        remove_block = [remove]*threads
        split_block = [args.split]*threads
        merge_block = [args.merge]*threads
        filter_block = [args.filter]*threads
        return zip(cns_blocks, select_block, remove_block, split_block, merge_block, filter_block, ass_block, ver_block)
    elif action == "coverage":        
        return zip(cns_blocks, samples_blocks, cols_block, ass_block, ver_block)
    elif action == "ploidy":
        return zip(cns_blocks, samples_blocks, cols_block, ass_block, ver_block)
    elif action == "bin":
        segs = load_regions(args.segments)
        segs_block = [segs]*threads
        fun_block = [args.aggregate]*threads
        return zip(cns_blocks, segs_block, fun_block, cols_block, ver_block)
    else:
        raise ValueError(f"Unknown action {action}")


def _process(action, cns_df, samples_df, cn_cols, assembly, args):   
    main_fun = _action_to_fun(action)
    blocks = abs(args.threads)
    samples_blocks = dataframe_array_split(samples_df, blocks)
    cns_blocks = [cns_df.query("sample_id in @block.index").reset_index(drop=True) for block in samples_blocks]
    zip_blocks = _get_blocks(action, cns_blocks, samples_blocks, cn_cols, assembly, blocks, args)
    if blocks == 1:
        return [main_fun(*list(*zip_blocks))]
    else:
        with Pool(blocks) as pool:
            log_info(args.verbose, f"Multiprocessing with {blocks} threads..")
            res_blocs = pool.starmap(main_fun, zip_blocks)        
            pool.close()   
            pool.join()    
        return res_blocs
    

def parse_cncols(cncols):
    if cncols != None:
        cncols = cncols.split(",")
        if len(cncols) == 1:
            cols_no = [cncols[0]]
        elif len(cncols) > 2:
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
    no_header = args.noheader
    no_sample = args.nosample
    subsplit = args.subsplit
    cncols = parse_cncols(args.cncols)

    # Read the input
    log_info(print_info, f"***** cns {action} *****")

    # Check the CNS file, may be omitted for segmentation without merging
    if action != "segment" or args.merge > 0:
        if not exists(cns_file_path):    
            raise ValueError(f"Copy number file {cns_file_path} not found.")
        else:
            log_info(print_info, f"CNS file at {cns_file_path}...")     
        
    cns_df = load_cns(cns_file_path, canonize=True, cn_columns=cncols, header=not no_header, no_sample=no_sample)
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
        res_dfs = _process(action, cns_df, samples_block, cn_columns, assembly, args)
        runtime = time.time() - start

        if print_info:
            print(f"Finished in {runtime:.3f} seconds. Writing to {out_file}...")        
            if args.time:  
                timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))  
                with open("./out/times.tsv", "a") as f:
                    f.write(f"{timestamp}\t{action}\t{args.threads}\t{out_file}\t{runtime}\n")

        for j in range(len(res_dfs)): 
            mode = "w" if i == 0 and j == 0 else "a"
            header = not no_header and i == 0 and j == 0
            res_df = res_dfs[j]
            if action == "segment":
                save_regions(res_df, out_file, change_coords=True, header=header, write_mode=mode)
            elif action in ["fill", "impute", "bin"]:
                save_cns(res_df, out_file, change_coords=True, no_sample=no_sample, header=header, write_mode=mode)
            else:
                if no_sample:
                    res_df.reset_index(drop=True)
                res_df.to_csv(out_file, sep="\t", index=True, header=header, mode=mode)

    log_info(print_info, "Done.")


if __name__ == "__main__":
    main()