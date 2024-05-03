#!/bin/bash

threads=8
preprocess_only=false

mkdir -p out
# TRACERx PCAWG TCGA_hg19 TCGA_hg38
for dataset in TRACERx PCAWG TCGA_hg19 TCGA_hg38;
do
    echo "Processing $dataset"
    ./data_preprocess.py ${dataset} --verbose
    if [ "$preprocess_only" = true ]; then
        continue      
    fi
    common_args="--threads $threads --verbose"
    cns fill "./out/${dataset}_cns_preprocess.tsv" --samples "./out/${dataset}_samples_preprocess.tsv" --out "./out/${dataset}_cns_fill.tsv" $common_args
    cns coverage "./out/${dataset}_cns_fill.tsv" --samples "./out/${dataset}_samples_preprocess.tsv" --out "./out/${dataset}_samples.tsv" $common_args
    cns impute "./out/${dataset}_cns_fill.tsv" --samples "./out/${dataset}_samples.tsv" --out "./out/${dataset}_cns_imp.tsv" $common_args
    cns ploidy "./out/${dataset}_cns_imp.tsv" --samples "./out/${dataset}_samples.tsv" --out "./out/${dataset}_samples.tsv" $common_args
done

