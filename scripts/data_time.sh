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
    common_args="--threads $threads --verbose"
    cns fill "${out}/${dataset}_cns_preprocess.tsv" --samples "${out}/${dataset}_samples_preprocess.tsv" --out "${temp_folder}/${dataset}_cns_fill.tsv" $common_args
    cns coverage "${out}/${dataset}_cns_fill.tsv" --samples "${out}/${dataset}_samples_preprocess.tsv" --out "${temp_folder}/${dataset}_samples.tsv" $common_args
    cns impute "${out}/${dataset}_cns_fill.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${temp_folder}/${dataset}_cns_imp.tsv" $common_args
    cns ploidy "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${temp_folder}/${dataset}_samples.tsv" $common_args
    cns bin "${out}/${dataset}_cns_imp.tsv" --out "${temp_folder}/${dataset}_bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000 $common_args
    cns bin "${out}/${dataset}_cns_imp.tsv" --select "${data}/COSMIC_consensus_genes.tsv" --out "${temp_folder}/${dataset}_bin_COSMIC.tsv" $common_args
done
rm -r $temp_folder
