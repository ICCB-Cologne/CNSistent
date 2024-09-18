#!/bin/bash

preprocess=true
process=true
postprocess=false

threads=30

data="../data"
out="../out"
cd "$(dirname "$0")" # Set path to the script's path

mkdir -p $out

# TRACERx_met TRACERx_prim PCAWG TCGA_hg19 TCGA_hg38
for dataset in TCGA_hg19;
do
    # Set assembly variable based on dataset
    if [ "$dataset" = "TCGA_hg38" ]; then
        assembly="hg38"
    else
        assembly="hg19"
    fi

    echo "Processing $dataset with assembly $assembly"    
    if [ "$preprocess" = true ]; then        
        echo "Preprocessing"
        ./data_preprocess.py ${dataset}
    fi
    if [ "$process" = true ]; then        
        common_args="--threads $threads --verbose --assembly $assembly"
        cns fill "${out}/${dataset}_cns_preprocess.tsv" --samples "${out}/${dataset}_samples_preprocess.tsv" --out "${out}/${dataset}_cns_fill.tsv" $common_args
        cns coverage "${out}/${dataset}_cns_fill.tsv" --samples "${out}/${dataset}_samples_preprocess.tsv" --out "${out}/${dataset}_samples.tsv" $common_args
        cns impute "${out}/${dataset}_cns_fill.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${out}/${dataset}_cns_imp.tsv" $common_args
        cns ploidy "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${out}/${dataset}_samples.tsv" $common_args
        cns signatures "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${out}/${dataset}_samples.tsv" $common_args
    fi
done

echo "Postprocessing"
if [ "$postprocess" = true ]; then
    ./data_postprocess.py
fi