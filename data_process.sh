#!/bin/bash

threads=8
preprocess_only=false

mkdir -p out
# TRACERx PCAWG TCGA_hg19 TCGA_hg38
for dataset in PCAWG;
do
    echo "Processing $dataset"
    ./data_preprocess.py ${dataset} --verbose
    if [ "$preprocess_only" = true ]; then
        continue      
    fi
    common_args="--threads $threads --verbose"
    ./cns.py fill "./out/${dataset}_cna_preprocess.tsv" --samples "./out/${dataset}_samples_preprocess.tsv" --out "./out/${dataset}_cna_fill.tsv" $common_args
    ./cns.py coverage "./out/${dataset}_cna_fill.tsv" --samples "./out/${dataset}_samples_preprocess.tsv" --out "./out/${dataset}_samples.tsv" $common_args
    ./cns.py impute "./out/${dataset}_cna_fill.tsv" --samples "./out/${dataset}_samples.tsv" --out "./out/${dataset}_cna_imp.tsv" $common_args
    ./cns.py ploidy "./out/${dataset}_cna_imp.tsv" --samples "./out/${dataset}_samples.tsv" --out "./out/${dataset}_samples.tsv" $common_args
done

