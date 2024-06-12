#!/bin/bash

threads=10
subsplit=10
data="../data"
out="../out"
cd "$(dirname "$0")" # Set path to the script's path

# TRACERx PCAWG TCGA_hg19 TCGA_hg38
for dataset in TCGA_hg19; 
do    
    echo "Processing $dataset"
    bin_shared_args="${out}/${dataset}_cns_imp.tsv --samples ${out}/${dataset}_samples.tsv  --verbose --threads $threads --subsplit $subsplit"
    for dist in 100000 300000 1000000; do
        let filter_dist=$dist/2
        cns cluster "${out}/${dataset}_cns_imp.tsv" --dist $dist --out "${out}/${dataset}_cluster_${dist}.tsv" --verbose
        cns bin $bin_shared_args --out "${out}/${dataset}_bin_cluster_${dist}.tsv" --select "${out}/${dataset}_cluster_${dist}.tsv" --remove gaps --filter $filter_dist
    done
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_10MB.tsv" --bins 10000000 --remove gaps --filter 1000000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_3MB.tsv" --bins 3000000 --remove gaps --filter 300000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_300KB.tsv" --bins 300000 --remove gaps --filter 30000
    cns bin $bin_shared_args --select "${data}/COSMIC_consensus_genes.tsv" --out "${out}/${dataset}_bin_COSMIC.tsv" --aggregate min
    cns bin $bin_shared_args --select "${data}/ENSEMBL_coding_genes.tsv" --out "${out}/${dataset}_bin_ENSEMBL.tsv" --aggregate min
    cns bin $bin_shared_args --select "arms" --out "${out}/${dataset}_bin_arms.tsv" --aggregate min --remove gaps --filter 1000000
    cns bin $bin_shared_args --select "bands" --out "${out}/${dataset}_bin_bands.tsv" --aggregate min --remove gaps --filter 100000
done
