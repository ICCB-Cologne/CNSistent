#!/bin/bash

data="../data"
out="../out"
cd "$(dirname "$0")" # Set path to the script's path

dataset=PCAWG
for threads in 1 2 4 8 16 32;
do
    echo "Processing $threads thread/s"
    common_args="--threads $threads --verbose"
    ./cns.py fill "${out}/${dataset}_cns_preprocess.tsv" --samples "${out}/${dataset}_samples_preprocess.tsv" --out "./temp/${dataset}_cns_fill.tsv" $common_args
    ./cns.py coverage "${out}/${dataset}_cns_fill.tsv" --samples "${out}/${dataset}_samples_preprocess.tsv" --out "./temp/${dataset}_samples.tsv" $common_args
    ./cns.py impute "${out}/${dataset}_cns_fill.tsv" --samples "${out}/${dataset}_samples.tsv" --out "./temp/${dataset}_cns_imp.tsv" $common_args
    ./cns.py ploidy "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "./temp/${dataset}_samples.tsv" $common_args
    ./cns.py bin "${out}/${dataset}_cns_imp.tsv" --out "./temp/${dataset}_bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000 $common_args
    ./cns.py bin "${out}/${dataset}_cns_imp.tsv" --select "${data}/COSMIC_consensus_genes.tsv" --out "./temp/${dataset}_bin_COSMIC.tsv" $common_args
done

