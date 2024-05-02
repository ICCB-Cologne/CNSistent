#!/bin/bash
preprocess=false

mkdir -p out
# TRACERx PCAWG TCGA_hg19 TCGA_hg38
dataset=PCAWG
for threads in 1 2 4 8 16 32;
do
    echo "Processing $threads thread/s"
    common_args="--threads $threads --verbose"
    ./cns.py fill "./out/${dataset}_cna_preprocess.tsv" --samples "./out/${dataset}_samples_preprocess.tsv" --out "./temp/${dataset}_cna_fill.tsv" $common_args
    ./cns.py coverage "./out/${dataset}_cna_fill.tsv" --samples "./out/${dataset}_samples_preprocess.tsv" --out "./temp/${dataset}_samples.tsv" $common_args
    ./cns.py impute "./out/${dataset}_cna_fill.tsv" --samples "./out/${dataset}_samples.tsv" --out "./temp/${dataset}_cna_imp.tsv" $common_args
    ./cns.py ploidy "./out/${dataset}_cna_imp.tsv" --samples "./out/${dataset}_samples.tsv" --out "./temp/${dataset}_samples.tsv" $common_args
    ./cns.py bin "./out/${dataset}_cna_imp.tsv" --out "./temp/${dataset}_bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000 $common_args
    ./cns.py bin "./out/${dataset}_cna_imp.tsv" --select "./data/COSMIC_consensus_genes.tsv" --out "./temp/${dataset}_bin_COSMIC.tsv" $common_args
done

