#!/bin/bash

# Function to run a command and compare its output
run_and_compare() {
  echo "Test: $1"

  # Run the command passed to the function
  eval "$1"

  # Extract the output file path
  # This splits the command by space and looks for the index of "--out" then gets the next item as the file path
  command_array=($1)
  for i in "${!command_array[@]}"; do
    if [[ "${command_array[$i]}" == "--out" ]]; then
      output_file="${command_array[$((i+1))]}"
      break
    fi
  done

  # Construct the path to the expected file
  expected_file="./out/$(basename "${output_file}")"

  # Compare the output with the expected output
  if ! diff "${output_file}" "${expected_file}" > /dev/null; then
    echo -e "\033[0;31m Test failed: Output of $1 does not match expected output ${expected_file}"
	diff "${output_file}" "${expected_file}"
  echo -e "Failed on command $1"
    exit 1
  fi
}

# Commands to run and test
commands=(
  "cns fill ./in/test_cns_double.tsv --samples ./in/test_samples.tsv --out ./temp/test_cns_fill.tsv --verbose"
  "cns fill ./in/test_cns_single.tsv --out ./temp/test_cns_single_fill.tsv --verbose"
  "cns impute ./temp/test_cns_fill.tsv --out ./temp/test_cns_imp.tsv --verbose"
  "cns impute ./temp/test_cns_single_fill.tsv --out ./temp/test_cns_single_imp.tsv --verbose"
  "cns coverage ./temp/test_cns_fill.tsv --out ./temp/test_sample_cover.tsv --verbose"
  "cns ploidy ./temp/test_cns_imp.tsv --samples ./in/test_samples.tsv --out ./temp/test_sample_ploidy.tsv --verbose"
  "cns breakage ./temp/test_cns_imp.tsv --samples ./in/test_samples.tsv --out ./temp/test_sample_breakage.tsv --verbose"
  "cns segment arms --out ./temp/test_segs_arms.bed --verbose"
  "cns segment bands --out ./temp/test_segs_bands.bed --verbose"
  "cns segment ./temp/test_cns_fill.tsv --split 1000000 --out ./temp/test_segs_1MB.bed --verbose"
  "cns segment ./temp/test_cns_fill.tsv --split 1000000 --out ./temp/test_segs_1MB_gaps.bed --remove gaps --filter 500000 --verbose"
  "cns segment arms --out ./temp/test_segs_arms_gaps.bed --remove gaps --filter 100000"
  "cns segment ../data/COSMIC_consensus_genes.bed --out ./temp/test_COSMIC_gaps.bed --remove gaps --filter 100000"
  "cns segment ./temp/test_cns_fill.tsv --merge 100000 --out ./temp/mcs_regions.bed --verbose "
  "cns aggregate ./temp/test_cns_fill.tsv --segments ./temp/test_segs_1MB.bed --out ./temp/test_cns_1MB.tsv --verbose"
  "cns aggregate ./temp/test_cns_imp.tsv --segments ./temp/test_COSMIC_gaps.bed --out ./temp/test_cns_COSMIC.tsv --verbose --how min"
  "cns aggregate ./in/test_cns_single.tsv --segments ./temp/test_segs_arms_gaps.bed --out ./temp/test_cns_arms.tsv --verbose"
)

# TODO: hg38 test

cd "$(dirname "$0")" # enter script dir
rm -r ./temp
mkdir ./temp
# Iterate over commands and run them
for cmd in "${commands[@]}"; do
  run_and_compare "$cmd"
done

echo "All tests passed successfully."
