#!/bin/bash

data="../data"
out="../out"
temp_folder="./temp"
cd "$(dirname "$0")" # Set path to the script's path

mkdir -p $temp_folder
dataset=PCAWG
for threads in 1 2 4 8 16 32;
do
    echo "Processing $threads thread/s"
    common_args="--threads $threads --verbose --time"
    cns align "${data}/${dataset}_cns_raw.tsv" --samples "${data}/${dataset}_samples_raw.tsv" --out "${temp_folder}/${dataset}_cns_align.tsv" $common_args
    cns coverage "${out}/${dataset}_cns_align.tsv" --samples "${data}/${dataset}_samples_raw.tsv" --out "${temp_folder}/${dataset}_samples.tsv" $common_args
    cns infer "${out}/${dataset}_cns_align.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${temp_folder}/${dataset}_cns_imp.tsv" $common_args
    cns ploidy "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${temp_folder}/${dataset}_samples.tsv" $common_args
    cns aggregate "${out}/${dataset}_cns_imp.tsv" --segments "${out}/segs_1MB.bed" --out "${temp_folder}/${dataset}_bin_1MB.tsv" $common_args
    cns aggregate "${out}/${dataset}_cns_imp.tsv" --segments "${out}/segs_COSMIC.bed" --out "${temp_folder}/${dataset}_bin_COSMIC.tsv" $common_args
done
rm -r $temp_folder
