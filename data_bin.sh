#!/bin/bash

threads=30

# TRACERx PCAWG TCGA_hg19 TCGA_hg38
for dataset in TCGA_hg19 TCGA_hg38 PCAWG TRACERx; 
do    
    echo "Processing $dataset"
    common_input="./out/${dataset}_cns_imp.tsv --samples ./out/${dataset}_samples.tsv"
    for dist in 100000 300000 1000000; do
        let filter_dist=$dist/2
        cns cluster $common_input --dist $dist --out "./out/${dataset}_cluster_${dist}.tsv" --verbose
        cns bin $common_input --out "./out/${dataset}_bin_cluster_${dist}.tsv" --select "./out/${dataset}_cluster_${dist}.tsv" --remove gaps --filter $filter_dist --threads $threads --verbose
    done
    cns bin $common_input --out "./out/${dataset}_bin_10MB.tsv" --bins 10000000 --remove gaps --filter 1000000 --threads $threads --verbose
    cns bin $common_input --out "./out/${dataset}_bin_3MB.tsv" --bins 3000000 --remove gaps --filter 300000 --threads $threads --verbose
    cns bin $common_input --out "./out/${dataset}_bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000 --threads $threads --verbose
    cns bin $common_input --out "./out/${dataset}_bin_300KB.tsv" --bins 300000 --remove gaps --filter 30000 --threads $threads --verbose
    cns bin $common_input --select "./data/COSMIC_consensus_genes.tsv" --out "./out/${dataset}_bin_COSMIC.tsv" --aggregate min --verbose --threads $threads
    cns bin $common_input --select "./data/ENSEMBL_coding_genes.tsv" --out "./out/${dataset}_bin_ENSEMBL.tsv" --aggregate min --verbose --threads $threads
    cns bin $common_input --select "arms" --out "./out/${dataset}_bin_arms.tsv" --aggregate min --verbose --threads $threads --remove gaps --filter 1000000
    cns bin $common_input --select "bands" --out "./out/${dataset}_bin_bands.tsv" --aggregate min --verbose --threads $threads --remove gaps --filter 100000
done
