#!/bin/bash

set -x

# set current directory to the location of the script
cd "$(dirname "$0")"

threads=30
subsplit=1
data="../data"
out="../out"
segment=true
cd "$(dirname "$0")" # Set path to the script's path

if [ $segment = true ]; then
    for dist in 1000000 500000 250000; do
        python ./data_cluster.py $dist
    done
    cns segment "dummy" --out "${out}/segs_10MB.bed" --split 10000000 --remove gaps --filter 1000000 
    cns segment "dummy" --out "${out}/segs_5MB.bed" --split 5000000 --remove gaps --filter 500000
    cns segment "dummy" --out "${out}/segs_3MB.bed" --split 3000000 --remove gaps --filter 300000  
    cns segment "dummy" --out "${out}/segs_2MB.bed" --split 2000000 --remove gaps --filter 200000
    cns segment "dummy" --out "${out}/segs_1MB.bed" --split 1000000 --remove gaps --filter 100000
    cns segment "dummy" --out "${out}/segs_500KB.bed" --split 500000 --remove gaps --filter 50000
    cns segment "dummy" --out "${out}/segs_250KB.bed" --split 250000 --remove gaps --filter 25000
    cns segment "dummy" --select "${data}/COSMIC_consensus_genes.bed" --out "${out}/segs_COSMIC.bed" 
    cns segment "dummy" --select "${data}/ENSEMBL_coding_genes.bed" --out "${out}/segs_ENSEMBL.bed"
    cns segment "dummy" --out "${out}/segs_whole.bed" --remove gaps --filter 1000000
    cns segment "dummy" --select "arms" --out "${out}/segs_arms.bed" --remove gaps --filter 1000000
    cns segment "dummy" --select "bands" --out "${out}/segs_bands.bed" --remove gaps --filter 100000
fi


# TRACERx PCAWG TCGA_hg19
for dataset in TRACERx PCAWG TCGA_hg19; 
do    
    echo "Processing $dataset"
    shared_args="${out}/${dataset}_cns_imp.tsv --samples ${out}/${dataset}_samples.tsv  --verbose --threads $threads --subsplit $subsplit"
    cns aggregate --segments "${out}/segs_whole.bed" --out "${out}/${dataset}_bin_whole.tsv" $shared_args
    cns aggregate --segments "${out}/segs_arms.bed" --out "${out}/${dataset}_bin_arms.tsv" $shared_args
    cns aggregate --segments "${out}/segs_bands.bed" --out "${out}/${dataset}_bin_bands.tsv" $shared_args
    cns aggregate --segments "${out}/segs_10MB.bed" --out "${out}/${dataset}_bin_10MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_5MB.bed" --out "${out}/${dataset}_bin_5MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_3MB.bed" --out "${out}/${dataset}_bin_3MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_2MB.bed" --out "${out}/${dataset}_bin_2MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_1MB.bed" --out "${out}/${dataset}_bin_1MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_500KB.bed" --out "${out}/${dataset}_bin_500KB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_250KB.bed" --out "${out}/${dataset}_bin_250KB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_COSMIC.bed" --out "${out}/${dataset}_bin_COSMIC.tsv" --how min $shared_args
    cns aggregate --segments "${out}/segs_ENSEMBL.bed" --out "${out}/${dataset}_bin_ENSEMBL.tsv" --how min $shared_args
    for dist in 1000000 500000 250000; do
        cns aggregate --segments "${out}/segs_merge_${dist}.bed" --out "${out}/${dataset}_bin_merge_${dist}.tsv" $shared_args 
    done
done

