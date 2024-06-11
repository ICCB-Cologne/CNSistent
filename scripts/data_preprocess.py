#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
from cns.utils.files import canonize_cns_df
from cns.utils.assemblies import get_assembly

_data_dir = "../data"
_out_dir = "../out"

def pcawg(print_debug=False):
    cns_raw_df = pd.read_csv(f"{_data_dir}/pcawg_cns_source.tsv", sep="\t")
    cns_raw_df = canonize_cns_df(cns_raw_df)
    if print_debug:
        print(cns_raw_df.head())
    specimen_df = pd.read_csv(f"{_data_dir}/pcawg_specimen_histology_August.tsv", sep="\t")
    sample_to_donor_df = specimen_df[["# icgc_specimen_id", "icgc_donor_id"]].rename(columns={"# icgc_specimen_id": "sample_id", "icgc_donor_id": "donor_id"})
    sample_to_donor_df = sample_to_donor_df.drop_duplicates().sort_values(by="sample_id").reset_index(drop=True)
    sample_to_donor_map = sample_to_donor_df.set_index("sample_id").to_dict()["donor_id"]
    if print_debug:
        print(sample_to_donor_df.head())

    donor_df = pd.read_csv(f"{_data_dir}/pcawg_donor_clinical_August2016.tsv", sep="\t")
    donor_to_sex_df = donor_df[["icgc_donor_id", "donor_sex"]].rename(columns={"icgc_donor_id": "donor_id"})
    donor_to_sex_df = donor_to_sex_df.drop_duplicates().sort_values(by="donor_id").reset_index(drop=True)
    donor_to_sex_map = donor_to_sex_df.set_index("donor_id").to_dict()["donor_sex"]
    if print_debug:
        print(donor_to_sex_df.head())
    ids = cns_raw_df["sample_id"].unique()
    donors = map(lambda x: sample_to_donor_map[x], ids)
    sexes = map(lambda x: donor_to_sex_map[x], donors)
    samples_df = pd.DataFrame([ids, sexes]).T
    samples_df.columns = ["sample_id", "sex"]
    samples_df = samples_df.set_index("sample_id")
    samples_df["sex"] = samples_df["sex"].replace({"male": "xy", "female": "xx"})
    # sanity check  
    assert len(cns_raw_df["sample_id"].unique()) == len(samples_df)

    # Add whitelist information
    whitelist = pd.read_csv(f"{_data_dir}/pcawg_supplementary_table_1.tsv", sep="\t")
    whitelist_sp = whitelist["icgc_specimen_id"].unique()
    assert(len(whitelist_sp) == 2583) # 2583 white-listed samples (https://www.nature.com/articles/s41586-020-1969-6)
    samples_df['whitelist'] = samples_df.index.isin(whitelist_sp)

    # Add type
    types = specimen_df[["# icgc_specimen_id", "histology_abbreviation"]].drop_duplicates().set_index("# icgc_specimen_id")
    types.index.name = "sample_id"
    samples_df["type"] = types

    # Add TCGA id
    ids = specimen_df[["# icgc_specimen_id", "submitted_donor_id"]].drop_duplicates().set_index("# icgc_specimen_id")
    ids.index.name = "sample_id"
    # set ids to empty string if does not start with TCGA
    ids["submitted_donor_id"] = np.where(ids["submitted_donor_id"].str.startswith("TCGA"), ids["submitted_donor_id"], "")
    samples_df["TCGA_id"] = ids

    # Add TCGA type
    msi_df = pd.read_csv(f"{_data_dir}/MS_analysis.PCAWG_release_v1.RIKEN.tsv", sep="\t")
    type_df = msi_df[["ID", "cancer"]].copy()
    TCGA_type_list  =[]
    for row in type_df.iterrows():
        id = row[1]["ID"]
        # find all the rows in specimen_df where id is in any column
        found = specimen_df.isin([id]).any(axis=1)
        if np.any(found):
            TCGA_type_list.append((specimen_df[found]["# icgc_specimen_id"].values[0], row[1]["cancer"]))# %%
    TCGA_type_df = pd.DataFrame(TCGA_type_list, columns=["sample_id", "cancer"]).set_index("sample_id")
    samples_df = samples_df.join(TCGA_type_df).rename(columns={"cancer": "TCGA_type"})

    if print_debug:
        print(samples_df.head())

    return cns_raw_df, samples_df


def tcga(hg_ver, print_debug=False):
    assembly = get_assembly(hg_ver)
    cns_source_df = pd.read_csv(f"{_data_dir}/tcga_{hg_ver}_cns_source.tsv", sep="\t")
    cns_raw_df = canonize_cns_df(cns_source_df, assembly=assembly)
    if print_debug:
        print(cns_raw_df.head())

    labels = f"summary.ascatv3TCGA.penalty70.{hg_ver}.tsv"
    specimen_df = pd.read_csv(f"{_data_dir}/{labels}", sep="\t")
    if print_debug:
        print(specimen_df.head())

    samples_df = specimen_df[["patient", "sex", "cancer_type"]].rename(columns={"patient": "sample_id", "cancer_type": "type"})
    samples_df["sex"] = samples_df["sex"].replace({"XY": "xy", "XX": "xx"})
    samples_df.set_index("sample_id", inplace=True)
    if print_debug:
        print(samples_df.head())

    return cns_raw_df, samples_df


def tracerx(print_debug=False):
    prim_cns_file = pd.read_csv(f"{_data_dir}/20220803_TxPri_mphase_by_sample_df.reduced.csv")
    if print_debug:
        print(prim_cns_file.head())

    cns_subset = prim_cns_file[["sample", "chr", "startpos", "endpos", "nMajor", "nMinor"]].copy()
    cns_subset["chr"] = "chr" + cns_subset["chr"].astype(str)
    cns_raw_df = canonize_cns_df(cns_subset, 2)
    if print_debug:
        print(cns_raw_df.head())

    labels = pd.read_csv(f"{_data_dir}/20221109_TRACERx421_all_patient_df.tsv", sep="\t")
    if print_debug:
        print(labels.head())

    left = prim_cns_file[["sample", "patient_id"]].drop_duplicates()
    right = labels[["cruk_id", "sex", "histology_multi_full"]].drop_duplicates()
    merged = pd.merge(left, right, left_on="patient_id", right_on="cruk_id")
    merged["sex"] = np.where(merged["sex"] == "Female", "xx", "xy")
    samples_df = merged[["sample", "sex", "histology_multi_full"]].rename(columns={"sample": "sample_id", "histology_multi_full": "type"})
    samples_df.set_index(["sample_id"], inplace=True)
    if print_debug:
        print(samples_df.head())

    return cns_raw_df, samples_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts source datasets into format usable by cns")
    parser.add_argument("dataset", type=str, help="Name of the dataset")  
    parser.add_argument("--verbose", action="store_true", help="Print debug information")
    args = parser.parse_args()
    dataset = args.dataset
    print_debug = args.verbose

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    if dataset == "PCAWG":
        cns_df, samples_df = pcawg(print_debug)
    elif dataset == "TCGA_hg19":
        cns_df, samples_df = tcga("hg19", print_debug)
    elif dataset == "TCGA_hg38":
        cns_df, samples_df = tcga("hg38", print_debug)
    elif dataset == "TRACERx":
        cns_df, samples_df = tracerx(print_debug)
    else:
        raise ValueError(f"Dataset {dataset} not recognized.")
    
    assert len(cns_df["sample_id"].unique()) == len(samples_df)
    samples_df.to_csv(f"{_out_dir}/{dataset}_samples_preprocess.tsv", sep="\t", index=True, header=True)
    cns_df.to_csv(f"{_out_dir}/{dataset}_cns_preprocess.tsv", sep="\t", index=False, header=True)