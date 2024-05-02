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
  expected_file="./tests/out/$(basename "${output_file}")"

  # Compare the output with the expected output
  if ! diff "${output_file}" "${expected_file}" > /dev/null; then
    echo -e "\033[0;31m Test failed: Output of $1 does not match expected output ${expected_file}"
	diff "${output_file}" "${expected_file}" 
    exit 1
  fi
}

# Commands to run and test
commands=(
  "python cns.py fill ./tests/in/test_cna_source.tsv --sample ./tests/in/test_sample_source.tsv --out ./temp/test_cna_fill.tsv"
  "python cns.py impute ./temp/test_cna_fill.tsv --out ./temp/test_cna_imp.tsv"
  "python cns.py coverage ./temp/test_cna_fill.tsv --out ./temp/test_sample_cover.tsv"
  "python cns.py ploidy ./tests/out/test_cna_imp.tsv --samples ./tests/in/test_sample_source.tsv --out ./tests/out/test_sample_ploidy.tsv"
  "python cns.py cluster ./temp/test_cna_fill.tsv --dist 100000 --out ./temp/mcs_regions.tsv"
  "python cns.py bin --select arms ./temp/test_cna_fill.tsv --out ./temp/test_cna_arms.tsv"
  "python cns.py bin --select bands ./temp/test_cna_fill.tsv --out ./temp/test_cna_bands.tsv"
  "python cns.py bin --bins 1000000 ./temp/test_cna_fill.tsv --out ./temp/test_cna_1MB.tsv"
  "python cns.py bin --bins 1000000 ./temp/test_cna_fill.tsv --out ./temp/test_cna_1MB_gaps.tsv --remove gaps --filter 500000"
  "python cns.py bin --select arms ./temp/test_cna_fill.tsv --out ./temp/test_segs_arms_gaps.tsv --remove gaps --filter 100000 --onlybins"
  "python cns.py bin ./temp/test_cna_fill.tsv --select ./temp/mcs_regions.tsv --out ./temp/test_cna_mcs.tsv"
)

rm -r ./temp
mkdir ./temp
# Iterate over commands and run them
for cmd in "${commands[@]}"; do
  run_and_compare "$cmd"
done

echo "All tests passed successfully."
