#!/bin/bash

threads=10
subsplit=5
data="../data"
out="../out"
segment=true
cd "$(dirname "$0")" # Set path to the script's path

if [ $segment = true ]; then
    for dist in 1000000 500000 250000; do
        python ./data_cluster.py $dist
    done
    cns segment "dummy" --out "${out}/segs_10MB.tsv" --split 10000000 --remove gaps --filter 1000000 
    cns segment "dummy" --out "${out}/segs_5MB.tsv" --split 5000000 --remove gaps --filter 500000
    cns segment "dummy" --out "${out}/segs_3MB.tsv" --split 3000000 --remove gaps --filter 300000  
    cns segment "dummy" --out "${out}/segs_2MB.tsv" --split 2000000 --remove gaps --filter 200000
    cns segment "dummy" --out "${out}/segs_1MB.tsv" --split 1000000 --remove gaps --filter 100000
    cns segment "dummy" --out "${out}/segs_500KB.tsv" --split 500000 --remove gaps --filter 50000
    cns segment "dummy" --out "${out}/segs_250KB.tsv" --split 250000 --remove gaps --filter 25000
    cns segment "dummy" --select "${data}/COSMIC_consensus_genes.tsv" --out "${out}/segs_COSMIC.tsv" 
    cns segment "dummy" --select "${data}/ENSEMBL_coding_genes.tsv" --out "${out}/segs_ENSEMBL.tsv"
    cns segment "dummy" --select "arms" --out "${out}/segs_arms.tsv" --remove gaps --filter 1000000
    cns segment "dummy" --select "bands" --out "${out}/segs_bands.tsv" --remove gaps --filter 100000
fi

# TRACERx PCAWG TCGA_hg19
for dataset in TRACERx PCAWG TCGA_hg19; 
do    
    echo "Processing $dataset"
    shared_args="${out}/${dataset}_cns_imp.tsv --samples ${out}/${dataset}_samples.tsv  --verbose --threads $threads --subsplit $subsplit"
    cns aggregate --segments "${out}/segs_arms.tsv" --out "${out}/${dataset}_bin_arms.tsv" $shared_args
    cns aggregate --segments "${out}/segs_bands.tsv" --out "${out}/${dataset}_bin_bands.tsv" $shared_args
    cns aggregate --segments "${out}/segs_10MB.tsv" --out "${out}/${dataset}_bin_10MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_5MB.tsv" --out "${out}/${dataset}_bin_5MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_3MB.tsv" --out "${out}/${dataset}_bin_3MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_2MB.tsv" --out "${out}/${dataset}_bin_2MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_1MB.tsv" --out "${out}/${dataset}_bin_1MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_500KB.tsv" --out "${out}/${dataset}_bin_500KB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_250KB.tsv" --out "${out}/${dataset}_bin_250KB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_COSMIC.tsv" --out "${out}/${dataset}_bin_COSMIC.tsv" --aggregate min $shared_args
    cns aggregate --segments "${out}/segs_ENSEMBL.tsv" --out "${out}/${dataset}_bin_ENSEMBL.tsv" --aggregate min $shared_args
    for dist in 1000000 500000 250000; do
        cns aggregate --segments "${out}/segs_merge_${dist}.tsv" --out "${out}/${dataset}_bin_merge_${dist}.tsv" $shared_args 
    done
done

