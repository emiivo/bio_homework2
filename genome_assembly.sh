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


# Create and activate ragtagEnv Conda environment, download RagTag

conda create -n ragtagEnv
conda activate ragtagEnv
#if RagTag is not installed, install it
if ! conda list ragtag | grep -q "ragtag"; then
  conda install -c bioconda ragtag
fi


# List of directories for SPAdes output
files_spades=(
  "$OUTPUT_DIR/ERR204044_trimmed_SPADES_OUT"
  "$OUTPUT_DIR/SRR15131330_trimmed_SPADES_OUT"
  "$OUTPUT_DIR/SRR18214264_trimmed_SPADES_OUT"
)

# RagTag for SPAdes
for i in "${!files_spades[@]}"; do
  fasta_file_spades="${files_spades[i]}/contigs.fasta"
  # Run ragtag.py scaffold with -u option
  ragtag.py scaffold ~/HW2/references/CP015498.fasta "$fasta_file_spades" -o "$OUTPUT_DIR/${files_names[i]}_spades_reordered"
done

files_abyss=(
  "$OUTPUT_DIR/ERR204044_alternative_contigs.fasta"
  "$OUTPUT_DIR/SRR15131330_alternative_contigs.fasta"
  "$OUTPUT_DIR/SRR18214264_alternative_contigs.fasta"
)

# Output file names
output_names_abyss_ragtag=(
  "ERR204044_abyss_reordered"
  "SRR15131330_abyss_reordered"
  "SRR18214264_abyss_reordered"
)

# RagTag for ABySS
for i in "${!files_abyss[@]}"; do
  ragtag.py scaffold ~/HW2/references/CP015498.fasta "${files_abyss[i]}" -o "$OUTPUT_DIR/${output_names_abyss_ragtag[i]}"
done

# Deactivate ragtagEnv environment
conda deactivate

# Using quast results select the best/better assembly for each sample.
# For ERR204044: the assembly made with ABySS, because it reaches a higher cumulative
# length and that means the assembly will be more complete and contiguous.

# For SRR15131330: the assembly made with SPAdes, because it reaches a higher cumulative
# length and that means the assembly will be more complete and contiguous. Also,
# assembly with ABySS almost has a linear relationship, which can indicate that the assembly consists 
# of numerous short contigs without substantial long-range connections or larger-scale structural information. 


# For SRR18214264: the assembly made with SPAdes, because it reaches a higher cumulative
# length much quicker and that means the assembly will be more complete.
 
