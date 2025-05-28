import numpy as np
import cns.data_utils as cdu
import subprocess
from os.path import join as pjoin
import os
import cns

# Get the patient ids in TRACERx
sample_df, _ = cdu.main_load("arms", dataset="TRACERx")
sample_ids = sample_df.index.unique()
# Patient ids are the first segment of sample_ids with the form CRUK0000_XXXXX
# The second segment refers to 
patient_ids = [sample_id.split("_")[0] for sample_id in sample_ids]
# Get the unique ids
patient_ids = sorted(list(set(patient_ids)))
####### TODO: Do you want all the patients or just one?
patient_ids = patient_ids[:1]

root_dir = os.path.abspath(os.path.join(cdu.out_path, "trees"))
os.makedirs(root_dir, exist_ok=True)

segmentations = ["arms", "20MB", "10MB", "5MB", "3MB", "2MB", "1MB", "500KB", "250KB", "100KB"]
for segs in segmentations:
    sample_df, cns_df = cdu.main_load(segs, dataset="TRACERx", use_filter=False)
    output_dir = os.path.join(root_dir, segs)
    os.makedirs(output_dir, exist_ok=True)
    for patient_id in patient_ids:
        # Output variables
        patient_dir = pjoin(output_dir, patient_id)
        os.makedirs(patient_dir, exist_ok=True)
        patient_dir = os.path.abspath(patient_dir)
        filename = pjoin(patient_dir, "cns.tsv")
        # Filter for the patient's samples and rename columns for MEDICC2
        patient_df = cns_df[cns_df["sample_id"].str.contains(patient_id)]
        patient_df = patient_df.rename(columns={"major_cn": "cn_a", "minor_cn":"cn_b"})
        patient_df.to_csv(filename, sep="\t", index=False)
        # Run medicc
        params = ["medicc2", filename, patient_dir, "-j", "10"]
        print(params)
        subprocess.Popen(params)