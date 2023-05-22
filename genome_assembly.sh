#!/bin/bash

INPUT_DIR=~/HW2/inputs
OUTPUT_DIR=~/HW2/outputs

# Define the input file pairs
input_files=(
  "$OUTPUT_DIR/SRR15131330_trimmed_2.fastq.gz $OUTPUT_DIR/SRR15131330_trimmed_1.fastq.gz"
  "$OUTPUT_DIR/ERR204044_trimmed_2.fastq.gz $OUTPUT_DIR/ERR204044_trimmed_1.fastq.gz"
  "$OUTPUT_DIR/SRR18214264_trimmed_2.fastq.gz $OUTPUT_DIR/SRR18214264_trimmed_1.fastq.gz"
)


# Loop through the input files
for files in "${input_files[@]}"; do
  # Split the file names into separate variables
  read -r input_file_1 input_file_2 <<<"$files"

  # Run SPAdes command
  spades.py --careful -o "${input_file_1%_*}_SPADES_OUT" -1 "$input_file_1" -2 "$input_file_2"
done

# Alternative assembly was created with ABySS through my own computer in browser, not by script.
# Alternative assembly can be found in ~/HW2/outputs

# Next, QUAST was used to evaluate assemblies made with SPAdes and ABySS.
# QUAST results can be found in ~/HW2/outputs/quast



