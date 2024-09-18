#!/bin/bash

cd "$(dirname "$0")" # Set path to the script's path

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
  "cns fill ./in/test_cns_source.tsv --sample ./in/test_sample_source.tsv --out ./temp/test_cns_fill.tsv"
  "cns impute ./temp/test_cns_fill.tsv --out ./temp/test_cns_imp.tsv"
  "cns coverage ./temp/test_cns_fill.tsv --out ./temp/test_sample_cover.tsv"
  "cns ploidy ./out/test_cns_imp.tsv --samples ./in/test_sample_source.tsv --out ./temp/test_sample_ploidy.tsv"
  "cns signatures ./out/test_cns_imp.tsv --samples ./in/test_sample_source.tsv --out ./temp/test_sample_signatures.tsv"
  "cns segment ./temp/test_cns_fill.tsv --merge 100000 --out ./temp/mcs_regions.tsv"
  "cns segment --select arms ./temp/test_cns_fill.tsv --out ./temp/test_segs_arms.tsv --threads 2"
  "cns segment --select bands ./temp/test_cns_fill.tsv --out ./temp/test_segs_bands.tsv --subsplit 2"
  "cns segment --split 1000000 ./temp/test_cns_fill.tsv --out ./temp/test_segs_1MB.tsv"
  "cns segment --split 1000000 ./temp/test_cns_fill.tsv --out ./temp/test_segs_1MB_gaps.tsv --remove gaps --filter 500000"
  "cns segment --select arms ./temp/test_cns_fill.tsv --out ./temp/test_segs_arms_gaps.tsv --remove gaps --filter 100000"
  "cns bin ./temp/test_cns_fill.tsv --segments ./temp/test_segs_1MB.tsv --out ./temp/test_cns_1MB.tsv"
)

# TODO: Test signle-columns BED file

rm -r ./temp
mkdir ./temp
# Iterate over commands and run them
for cmd in "${commands[@]}"; do
  run_and_compare "$cmd"
done

echo "All tests passed successfully."
