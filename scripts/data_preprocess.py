#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
from cns.utils.canonization import canonize_cns_df
from cns.utils.assemblies import get_assembly
from cns.data_utils import load_cns_out, load_samples_out, out_path, data_path
from cns.utils.files import save_cns, save_samples
from cns.utils.logging import log_info


def pcawg(print_info=False):
    cns_raw_df = pd.read_csv(f"{data_path}/pcawg_cns_source.tsv", sep="\t")
    cns_raw_df = canonize_cns_df(cns_raw_df)
    specimen_df = pd.read_csv(f"{data_path}/pcawg_specimen_histology_August.tsv", sep="\t")
    rename_map = {"# icgc_specimen_id": "sample_id", "icgc_donor_id": "donor_id"}
    sample_to_donor_df = specimen_df[["# icgc_specimen_id", "icgc_donor_id"]].rename(columns=rename_map)
    sample_to_donor_df = sample_to_donor_df.drop_duplicates().sort_values(by="sample_id").reset_index(drop=True)
    sample_to_donor_map = sample_to_donor_df.set_index("sample_id").to_dict()["donor_id"]

    donor_df = pd.read_csv(f"{data_path}/pcawg_donor_clinical_August2016.tsv", sep="\t")
    donor_to_sex_df = donor_df[["icgc_donor_id", "donor_sex"]].rename(columns={"icgc_donor_id": "donor_id"})
    donor_to_sex_df = donor_to_sex_df.drop_duplicates().sort_values(by="donor_id").reset_index(drop=True)
    donor_to_sex_map = donor_to_sex_df.set_index("donor_id").to_dict()["donor_sex"]

    ids = cns_raw_df["sample_id"].unique()
    donors = map(lambda x: sample_to_donor_map[x], ids)
    sexes = map(lambda x: donor_to_sex_map[x], donors)
    samples_df = pd.DataFrame([ids, sexes]).T
    samples_df.columns = ["sample_id", "sex"]
    samples_df = samples_df.set_index("sample_id")
    samples_df["sex"] = samples_df["sex"].replace({"male": "xy", "female": "xx"})
    # Correct for these 6 as they have calls for chrY
    samples_df.loc[["SP107470", "SP107451", "SP107557", "SP107449", "SP107448", "SP107446"]] = "xy"
    # sanity check
    assert len(cns_raw_df["sample_id"].unique()) == len(samples_df)

    # Add whitelist information
    whitelist = pd.read_csv(f"{data_path}/pcawg_supplementary_table_1.tsv", sep="\t")
    whitelist_sp = whitelist["icgc_specimen_id"].unique()
    assert len(whitelist_sp) == 2583  # 2583 white-listed samples (https://www.nature.com/articles/s41586-020-1969-6)
    samples_df["whitelist"] = samples_df.index.isin(whitelist_sp)

    # Add type
    types = (
        specimen_df[["# icgc_specimen_id", "histology_abbreviation"]].drop_duplicates().set_index("# icgc_specimen_id")
    )
    hist2 = (
        specimen_df[["# icgc_specimen_id", "histology_tier2"]].drop_duplicates().set_index("# icgc_specimen_id")
    )
    types.index.name = "sample_id"
    hist2.index.name = "sample_id"
    samples_df["type"] = types
    samples_df["histology_tier2"] = hist2

    # Add TCGA id
    ids = specimen_df[["# icgc_specimen_id", "submitted_donor_id"]].drop_duplicates().set_index("# icgc_specimen_id")
    ids.index.name = "sample_id"
    # set ids to empty string if does not start with TCGA
    ids["submitted_donor_id"] = np.where(
        ids["submitted_donor_id"].str.startswith("TCGA"), ids["submitted_donor_id"], ""
    )
    samples_df["TCGA_id"] = ids

    # Add TCGA type
    msi_df = pd.read_csv(f"{data_path}/MS_analysis.PCAWG_release_v1.RIKEN.tsv", sep="\t")
    type_df = msi_df[["ID", "cancer"]].copy()
    TCGA_type_list = []
    for row in type_df.iterrows():
        id = row[1]["ID"]
        # find all the rows in specimen_df where id is in any column
        found = specimen_df.isin([id]).any(axis=1)
        if np.any(found):
            TCGA_type_list.append((specimen_df[found]["# icgc_specimen_id"].values[0], row[1]["cancer"]))  # %%
    TCGA_type_df = pd.DataFrame(TCGA_type_list, columns=["sample_id", "cancer"]).set_index("sample_id")
    samples_df = samples_df.join(TCGA_type_df).rename(columns={"cancer": "TCGA_type"})
    samples_df.loc["SP116919", "TCGA_type"] = "CMDI" # Not matched, correct manually

    return cns_raw_df, samples_df


def tcga(hg_ver, print_info=False):
    assembly = get_assembly(hg_ver)
    cns_source_df = pd.read_csv(f"{data_path}/tcga_{hg_ver}_cns_source.tsv", sep="\t")
    cns_raw_df = canonize_cns_df(cns_source_df, assembly=assembly)

    labels = f"summary.ascatv3TCGA.penalty70.{hg_ver}.tsv"
    specimen_df = pd.read_csv(f"{data_path}/{labels}", sep="\t")

    samples_df = specimen_df[["name", "sex", "cancer_type"]].rename(
        columns={"name": "sample_id", "cancer_type": "type"}
    )
    samples_df["sex"] = samples_df["sex"].replace({"XY": "xy", "XX": "xx"})
    samples_df.set_index("sample_id", inplace=True)

    return cns_raw_df, samples_df


def tracerx(primary, print_info=False):
    if primary:
        cns_file = pd.read_csv(f"{data_path}/20220803_TxPri_mphase_by_sample_df.reduced.csv")
    else:
        cns_file = pd.read_csv(f"{data_path}/20220807_TxMets_prim_and_met_mphase_by_sample_df.reduced.csv")

    cns_subset = cns_file[["sample", "chr", "startpos", "endpos", "nMajor", "nMinor"]].copy()
    cns_subset["chr"] = "chr" + cns_subset["chr"].astype(str)
    cns_raw_df = canonize_cns_df(cns_subset)

    labels = pd.read_csv(f"{data_path}/20221109_TRACERx421_all_patient_df.tsv", sep="\t")

    left = cns_file[["sample", "patient_id"]].drop_duplicates()
    right = labels[["cruk_id", "sex", "histology_multi_full"]].drop_duplicates()
    merged = pd.merge(left, right, left_on="patient_id", right_on="cruk_id")
    merged["sex"] = np.where(merged["sex"] == "Female", "xx", "xy")
    samples_df = merged[["sample", "sex", "histology_multi_full"]].rename(
        columns={"sample": "sample_id", "histology_multi_full": "type"}
    )
    samples_df.set_index(["sample_id"], inplace=True)

    return cns_raw_df, samples_df


def get_processed_data(dataset, print_info=False):
    if dataset == "PCAWG":
        return pcawg(print_info)
    elif dataset == "TCGA_hg19":
        return tcga("hg19", print_info)
    elif dataset == "TCGA_hg38":
        return tcga("hg38", print_info)
    elif dataset == "TRACERx_prim":
        return tracerx(True, print_info)
    elif dataset == "TRACERx_met":
        return tracerx(False, print_info)
    else:
        raise ValueError(f"Dataset {dataset} not recognized.")


def merge_TRACERx_samples(print_info=False):
    # Load the sample data
    prim_samples = load_samples_out(f"TRACERx_prim_samples_preprocess.tsv")
    met_samples = load_samples_out(f"TRACERx_met_samples_preprocess.tsv")

    # Merge the DataFrames on the common key 'sample_id'
    common_samples = pd.merge(prim_samples, met_samples, on="sample_id", suffixes=("_prim", "_met")).index.unique()

    # Filtered indices
    log_info(print_info, f"TRACERx primary samples: {len(prim_samples)}")
    log_info(print_info, f"TRACERx metastatic samples: {len(met_samples)}")
    log_info(print_info, f"TRACERx all samples: {len(prim_samples) + len(met_samples)}")
    log_info(print_info, f"TRACERx common samples: {len(common_samples)}")
    met_only = met_samples[~met_samples.index.isin(common_samples)].copy()

    # dataset label
    prim_samples["TRACERx_set"] = "primary"
    prim_samples[prim_samples.isin(common_samples)]["TRACERx_set"] = "both"
    met_only["TRACERx_set"] = "metastatic"

    all_df = pd.concat([prim_samples, met_only], axis=0)
    log_info(print_info, f"TRACERx all samples after merge: {len(all_df)}")

    # Save the merged DataFrame
    save_samples(all_df, f"{out_path}/TRACERx_samples_preprocess.tsv")


def merge_TRACERx_cns(print_info=False, filled=False):
    prim_cns = load_cns_out(f"TRACERx_prim_cns_preprocess.tsv", raw=True)
    prim_cns.set_index(["sample_id"], inplace=True)
    met_cns = load_cns_out(f"TRACERx_met_cns_preprocess.tsv", raw=True)
    met_cns.set_index(["sample_id"], inplace=True)

    common_samples = list(set(prim_cns.index).intersection(set(met_cns.index)))
    filtered_met_cns = met_cns[~met_cns.index.isin(common_samples)].reset_index()
    prim_cns.reset_index(inplace=True)

    res_df = pd.concat([prim_cns, filtered_met_cns]).sort_values(["sample_id", "chrom", "start"])
    save_cns(res_df, f"{out_path}/TRACERx_cns_preprocess.tsv")


def remove_PCAWG_from_TCGA(print_info=False):
    # Load the sample data
    TCGA_samples = load_samples_out("TCGA_hg19_samples_preprocess.tsv")
    PCAWG_samples = load_samples_out("PCAWG_samples_preprocess.tsv")

    overlap_with_tcga = PCAWG_samples["TCGA_id"].dropna().unique()
    log_info(print_info, f"Total TCGA samples: {len(TCGA_samples)}")
    log_info(print_info, f"Overlapping samples with PCAWG from total: {len(overlap_with_tcga)}")
    TCGA_samples = TCGA_samples.query("sample_id not in @overlap_with_tcga")
    log_info(print_info, f"Total TCGA samples after filtering: {len(TCGA_samples)}")

    # Save the merged DataFrame
    save_samples(TCGA_samples, f"{out_path}/TCGA_hg19_samples_preprocess.tsv")


def filter_TCGA_CNS(print_info=False):
    # Load the sample data
    TCGA_cns = load_cns_out("TCGA_hg19_cns_preprocess.tsv", raw=True)
    TCGA_samples = load_samples_out("TCGA_hg19_samples_preprocess.tsv")

    # Filter the CNS data
    TCGA_cns = TCGA_cns[TCGA_cns["sample_id"].isin(TCGA_samples.index)]

    # Save the filtered CNS data
    save_cns(TCGA_cns, f"{out_path}/TCGA_hg19_cns_preprocess.tsv")


def remove_not_whitelist_PCAWG(print_info=False):
    # Load the sample data
    PCAWG_samples = load_samples_out("PCAWG_samples_preprocess.tsv")

    blacklist_count = len(PCAWG_samples) - PCAWG_samples["whitelist"].sum()
    log_info(print_info, f"Total PCAWG samples: {len(PCAWG_samples)}")
    log_info(print_info, f"Blacklisted samples: {blacklist_count}")

    # Filter the samples
    PCAWG_samples = PCAWG_samples[PCAWG_samples["whitelist"]]
    PCAWG_samples.drop(columns=["whitelist"], inplace=True)
    log_info(print_info, f"Total PCAWG samples after filtering: {len(PCAWG_samples)}")

    # Save the filtered samples
    save_samples(PCAWG_samples, f"{out_path}/PCAWG_samples_preprocess.tsv")


def filter_PCAWG_CNS(print_info=False):
    # Load the CNS data
    PCAWG_cns = load_cns_out("PCAWG_cns_preprocess.tsv")
    PCAWG_samples = load_samples_out("PCAWG_samples_preprocess.tsv")

    # Filter the CNS data
    PCAWG_cns = PCAWG_cns[PCAWG_cns["sample_id"].isin(PCAWG_samples.index)]

    # Save the filtered CNS data
    save_cns(PCAWG_cns, f"{out_path}/PCAWG_cns_preprocess.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts source datasets into format usable by cns")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print debug information")
    args = parser.parse_args()
    should_print = args.verbose
    log_fun = lambda x: log_info(should_print, x)

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for dataset in ["PCAWG", "TCGA_hg19", "TCGA_hg38", "TRACERx_prim", "TRACERx_met"]:
        log_fun(f"Preprocessing {dataset} data")
        cns_df, samples_df = get_processed_data(dataset, should_print)
        assert len(cns_df["sample_id"].unique()) == len(samples_df)
        samples_df.to_csv(f"{out_path}/{dataset}_samples_preprocess.tsv", sep="\t", index=True, header=True)
        cns_df.to_csv(f"{out_path}/{dataset}_cns_preprocess.tsv", sep="\t", index=False, header=True)

    log_fun("Merging samples for TRACERx...")
    merge_TRACERx_samples(should_print)
    log_fun("Merging filled CNS for TRACERx...")
    merge_TRACERx_cns(should_print, True)
    log_fun("Merging imputed CNS for TRACERx...")
    merge_TRACERx_cns(should_print, False)
    log_fun("Filtering TCGA samples by PCAWG samples...")
    remove_PCAWG_from_TCGA(should_print)
    filter_TCGA_CNS(should_print)
    log_fun("Filtering PCAWG samples by whitelist...")
    remove_not_whitelist_PCAWG(should_print)
    filter_PCAWG_CNS(should_print)
    log_fun("Done!")
