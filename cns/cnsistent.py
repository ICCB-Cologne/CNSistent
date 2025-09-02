#!/usr/bin/env python

from multiprocessing import Pool
import time
import argparse
from os.path import exists

from cns.utils.files import obtain_segments
from cns.utils.selection import dataframe_array_split
from cns.utils.files import *
from cns.pipelines import *
from cns.utils.selection import dataframe_array_split


def _add_sp_args(action, parser):
    parser.add_argument("data", type=str, help="path to a TSV file with copy number segments, potentially with multiple samples.")
    parser.add_argument("--samples", type=str, help="path to the samples file", required=False, default="")
    parser.add_argument("--out", type=str, help="output file path", required=False, default="./cns.out.tsv")
    parser.add_argument("--assembly", type=str, help="assembly to use. One of: hg19, hg38.", required=False, default="hg19")
    parser.add_argument(
        "--cncols",
        type=str,
        help="The name of either a single CN column or two comma separated columns. E.g. 'total_cn' or 'major_cn,minor_cn'.",
        required=False,
        default=None,
    )
    parser.add_argument("--threads", type=int, help="number of threads to use", required=False, default=1)
    parser.add_argument(
        "--subsplit",
        type=int,
        help="will split the processing into chunks to lower memory needs",
        required=False,
        default=1
    )
    parser.add_argument("--verbose", help="print progress to console", action="store_true")
    parser.add_argument("--time", help="save runtime info", action="store_true")

    if action in ["infer", "impute"]:
        parser.add_argument(
            "--method",
            type=str,
            help='Inference method to use. Options are "extend", "diploid", or "zero". Default is "extend".',
            required=False,
            default="extend"
        )

    if action in ["align", "impute"]:
        parser.add_argument(
            "--add-missing-chroms",
            type=bool,
            help="If True, adds missing chromosomes to the data. Default is True.",
            required=False,
            default=True
        )

    if action in ["coverage", "ploidy", "breakage", "aggregate"]:
        parser.add_argument(
            "--segments",
            type=str,
            help="Either a path to a segmentation file or a predefined segment type (e.g., 'whole', 'arms', 'bands').",
            required=False,
            default="whole",
        )

    if action in ["segment"]:
        parser.add_argument(
            "--split",
            type=int,
            help="Distance in which regions should be, can be a positive integer or negative for no splitting (whole regions).",
            required=False,
            default=-1,
        )
        parser.add_argument(
            "--remove",
            type=str,
            help="Removed the regions after selection, before segmentation, can be either 'gaps', path to a BED file, or empty.",
            required=False,
            default="",
        )
        parser.add_argument(
            "--filter",
            type=int,
            help="If set, regions smaller than the given size are excluded from selection and gaps. If negative, no filtering is done.",
            required=False,
            default=-1,
        )
        parser.add_argument(
            "--merge",
            type=int,
            help="Maximum distance between breakpoint clusters for breakpoint merging. If negative, no breakpoints are merged.",
            required=False,
            default=-1,
        )

        
    if action == "aggregate":
        parser.add_argument(
            "--how",
            type=str,
            help="The aggregation function, one of ['min', 'max', 'mean']",
            required=False,
            default="mean",
        )


def _get_version():
    try:
        with open("pyproject.toml", "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("version"):
                    # Handles: version = "1.2.3"
                    return line.split("=")[1].strip().strip('"').strip("'")
    except Exception:
        return "unknown"


def _parse_args():
    # Top-level parser
    parser = argparse.ArgumentParser(description="Impute missing values in CNS data")
    subparsers = parser.add_subparsers(dest="action", help="cns action to perform. One of:")
    # Parse version from pyproject.toml

    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {_get_version()}")

    sp_dict = {}
    sp_dict["align"] = subparsers.add_parser("align", help=f"Adds Nan regions to the CNS data to match the assembly.")
    sp_dict["infer"] = subparsers.add_parser("infer", help=f"Infers values for NaNs in the CNS data.")
    sp_dict["impute"] = subparsers.add_parser("impute", help=f"Imputes missing values in the CNS data. (combines align and infer)")
    sp_dict["coverage"] = subparsers.add_parser("coverage", help=f"Calculates coverage for aligned (but not imputed) CNS data." )
    sp_dict["ploidy"] = subparsers.add_parser("ploidy", help=f"Conducts breakpoint analysis for CNS data (NaNs are ignored).")
    sp_dict["breakage"] = subparsers.add_parser("breakage", help=f"Extracts basal CN signatures from CNS data (NaNs are ignored).")
    sp_dict["segment"] = subparsers.add_parser("segment", help=f"Calculates segmentation regions for CNS data.")
    sp_dict["aggregate"] = subparsers.add_parser("aggregate", help=f"Aggregate copy numbers across segments to match provided segments.")
    for action, sp in sp_dict.items():
        _add_sp_args(action=action, parser=sp)

    args = parser.parse_args()
    if args.action is None:
        parser.print_help()
        exit(1)
    if args.action not in sp_dict:
        raise ValueError(f"Action {args.action} not recognized.")

    if args.action == "segment":
        if args.threads > 1:
            print("segmentation is not data parallelizable, --threads option will be ignored.")
            args.threads = 1
        if args.subsplit > 1:
            print("segmentation is not data parallelizable, --subsplit option will be ignored.")
            args.subsplit = 1
    else:
        if not exists(args.data):
           raise ValueError(f"Data file {args.data} not found.")

    if args.threads <= 0:
        raise ValueError("The --threads option must be greater than 0.")

    if args.subsplit <= 0:
        raise ValueError("The --subsplit option must be greater than 0.")

    return args


def _action_to_fun(action):
    if action == "align":
        return main_align
    elif action == "infer":
        return main_infer
    elif action == "impute":
        return main_impute
    elif action == "coverage":
        return main_coverage
    elif action == "ploidy":
        return main_ploidy
    elif action == "breakage":
        return main_breakage
    elif action == "aggregate":
        return main_aggregate
    else:
        raise ValueError(f"Action {action} not recognized.")


def _get_blocks(action, input_block, samples_blocks, cn_cols, segs_block, assembly, args):
    block_count = len(input_block)
    # Apply process_block to each pair of blocks
    ass_block = [assembly] * block_count
    ver_block = [False] * block_count
    ver_block[-1] = args.verbose
    cols_block = [cn_cols] * block_count
    if action == "infer":        
        method_block = [args.method] * block_count
        return zip(input_block, samples_blocks, method_block, cols_block, ver_block)
    if action == "align":
        add_missing = [args.add_missing_chroms] * block_count
        return zip(input_block, samples_blocks, cols_block, add_missing, ass_block, ver_block)
    if action == "impute":
        method_block = [args.method] * block_count
        add_missing = [args.add_missing_chroms] * block_count
        return zip(input_block, samples_blocks, cols_block, method_block, add_missing, ass_block, ver_block)
    elif action in ["coverage", "ploidy", "breakage"]:
        if segs_block is None:
            raise ValueError("Segmentation blocks must be provided for this action.")
        return zip(input_block, samples_blocks, cols_block, segs_block, ass_block, ver_block)
    elif action == "aggregate":
        if segs_block is None:
            raise ValueError("Segmentation blocks must be provided for this action.")
        fun_block = [args.how] * block_count
        return zip(input_block, segs_block, fun_block, cols_block, ver_block)
    else:
        raise ValueError(f"Unknown action {action}")


def _process(action, cns_df, samples_df, cn_cols, select_segs, assembly, args):
    main_fun = _action_to_fun(action)
    threads = abs(args.threads)
    samples_blocks = dataframe_array_split(samples_df, threads)
    cns_blocks = [cns_df.query("sample_id in @block.index").reset_index(drop=True) for block in samples_blocks]
    segs_blocks = [select_segs] * len(cns_blocks) if select_segs is not None else None
    zip_blocks = _get_blocks(action, cns_blocks, samples_blocks, cn_cols, segs_blocks, assembly, args)
    if threads == 1:
        return [main_fun(*list(*zip_blocks))]
    else:
        with Pool(threads) as pool:
            log_info(args.verbose, f"Multiprocessing with {threads} threads..")
            res_blocs = pool.starmap(main_fun, zip_blocks)
            pool.close()
            pool.join()
        return res_blocs


def _parse_cncols(cncols):
    if cncols != None:
        cncols = cncols.split(",")
        if len(cncols) > 2:
            raise ValueError("Only one or two columns can be specified.")
    return cncols


def _save_time(action, out_file, runtime, start, threads):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start))
        filepath = "../out/times.tsv"
        write_mode = "a" if exists(filepath) else "w"
        with open(filepath, write_mode) as f:
            f.write(f"{timestamp}\t{action}\t{threads}\t{out_file}\t{runtime}\n")


def main():
    args = _parse_args()
    action = args.action
    assembly = get_assembly(args.assembly)
    print_info = args.verbose
    in_cols = _parse_cncols(args.cncols)
    out_file = args.out

    log_info(print_info, f"***** cns {action} *****")

    # Process segments
    if action in ["segment"]:
        input_segs = obtain_segments(args.data, in_cols, assembly, print_info)        
        remove_regs = obtain_segments(args.remove, in_cols, assembly, print_info)   
        res_segs = main_segment(input_segs, remove_regs, args.split, args.merge, args.filter, print_info)
        save_segments(res_segs, out_file)
        return

    log_info(print_info, f"Loading CNS input file {args.data}...")
    input_data = load_cns(args.data, cn_columns=in_cols, assembly=assembly, print_info=print_info)
    cn_columns = get_cn_cols(input_data, in_cols)
    if args.samples == "":
        samples_df = samples_df_from_cns_df(input_data, False)
    else:
        samples_df = load_samples(args.samples)
    samples_df = fill_sex_if_missing(input_data, samples_df)
    samples_blocks = dataframe_array_split(samples_df, args.subsplit)
    select_segs = obtain_segments(args.segments, in_cols, assembly, print_info) if "segments" in args else None

    # Process blocks
    for i in range(args.subsplit):
        log_info(print_info, f"Processing block {i+1}/{args.subsplit}...")

        samples_block = samples_blocks[i]

        # Perform the action
        start = time.time()
        
        res_list = _process(action, input_data, samples_block, cn_columns, select_segs, assembly, args)
        runtime = time.time() - start

        if print_info:
            print(f"Finished in {runtime:.3f} seconds. Writing to {out_file} ...")
            if args.time:
                _save_time(action, out_file, runtime, start, args.threads)

        for j in range(len(res_list)):
            mode = "w" if i == 0 and j == 0 else "a"
            res = res_list[j]
            if action in ["align", "infer", "impute", "aggregate"]:
                save_cns(res, out_file, change_coords=True, mode=mode)
            elif action in ["coverage", "ploidy", "breakage"]:
                save_samples(res, out_file, mode=mode)
            else:
                raise ValueError(f"Unknown action {action}")

    log_info(print_info, "Done.")


if __name__ == "__main__":
    main()
