#!/bin/bash

threads=10
subsplit=5
data="../data"
out="../out"
cluster=true
cd "$(dirname "$0")" # Set path to the script's path

if [ $cluster = true ]; then
    for dist in 1000000 500000 250000; do
        python ./data_cluster.py $dist
    done   
fi

# TRACERx PCAWG TCGA_hg19 TCGA_hg38
for dataset in TRACERx PCAWG TCGA_hg19 TCGA_hg38; 
do    
    echo "Processing $dataset"
    bin_shared_args="${out}/${dataset}_cns_imp.tsv --samples ${out}/${dataset}_samples.tsv  --verbose --threads $threads --subsplit $subsplit"
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_10MB.tsv" --bins 10000000 --remove gaps --filter 1000000 
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_5MB.tsv" --bins 5000000 --remove gaps --filter 500000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_3MB.tsv" --bins 3000000 --remove gaps --filter 300000  
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_2MB.tsv" --bins 2000000 --remove gaps --filter 200000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_1MB.tsv" --bins 1000000 --remove gaps --filter 100000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_500KB.tsv" --bins 500000 --remove gaps --filter 50000
    cns bin $bin_shared_args --out "${out}/${dataset}_bin_250KB.tsv" --bins 250000 --remove gaps --filter 25000
    cns bin $bin_shared_args --select "${data}/COSMIC_consensus_genes.tsv" --out "${out}/${dataset}_bin_COSMIC.tsv" --aggregate min
    cns bin $bin_shared_args --select "${data}/ENSEMBL_coding_genes.tsv" --out "${out}/${dataset}_bin_ENSEMBL.tsv" --aggregate min
    cns bin $bin_shared_args --select "arms" --out "${out}/${dataset}_bin_arms.tsv" --aggregate min --remove gaps --filter 1000000
    cns bin $bin_shared_args --select "bands" --out "${out}/${dataset}_bin_bands.tsv" --aggregate min --remove gaps --filter 100000
    for dist in 1000000 500000 250000; do
        let filter_dist=$dist/10
        cns bin $bin_shared_args --out "${out}/${dataset}_bin_cluster_${dist}.tsv" --select "${out}/joint_clust_${dist}.tsv" --remove gaps --filter $filter_dist
    done
done

