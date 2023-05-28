#!/bin/bash

INPUT_DIR=~/HW2/inputs
OUTPUT_DIR=~/HW2/outputs
RESULTS_DIR=~/HW2/results

# -----------------------------------
# First part of the script is genome mapping with abyss.
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

# -----------------------------------------
# Alternative assembly was created with ABySS through my own computer in browser, not by script.
# Alternative assembly can be found in ~/HW2/outputs


# ------------------------------------------
# Next, QUAST was used to evaluate assemblies made with SPAdes and ABySS.
# QUAST results can be found in ~/HW2/outputs/quast


# ------------------------------------------
# Next part of the script is making scaffolds (orienating contigs) using ragtag. 
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

# ----------------------------------
# Using quast results select the best/better assembly for each sample:
# For ERR204044: the assembly made with ABySS, because it reaches a higher cumulative
# length and that means the assembly will be more complete and contiguous.

# For SRR15131330: the assembly made with SPAdes, because it reaches a higher cumulative
# length and that means the assembly will be more complete and contiguous. Also,
# assembly with ABySS almost has a linear relationship, which can indicate that the assembly consists 
# of numerous short contigs without substantial long-range connections or larger-scale structural information. 


# For SRR18214264: the assembly made with SPAdes, because it reaches a higher cumulative
# length much quicker and that means the assembly will be more complete.

# ---------------------------------------------
# Use BWA for mapping. I decided to use BWA because it supports the mapping of
# reads to a single reference genome without the need for transcriptome-based 
# mapping or splicing-aware alignment.


# Define output directory
OUTPUT_DIR=~/HW2/outputs
INPUT_DIR=~/HW2/inputs

# Define abyss assembly file names
abyss_assembly_files=(
  "ERR204044_alternative_contigs.fasta"
  "SRR15131330_alternative_contigs.fasta"
  "SRR18214264_alternative_contigs.fasta"
)

# Define read file names
read_files=(
  "ERR204044"
  "SRR15131330"
  "SRR18214264"
)

# Map each read file to its corresponding alternative assembly file
for ((i=0; i<${#abyss_assembly_files[@]}; i++)); do
  assembly_file="${abyss_assembly_files[i]}"
  read_file="${read_files[i]}"
  sam_output="$OUTPUT_DIR/${assembly_file%.fasta}_$read_file.abyss.sam"
  bwa mem -t 6 "$OUTPUT_DIR/$assembly_file" "$INPUT_DIR/$read_file"_1.fastq.gz "$INPUT_DIR/$read_file"_2.fastq.gz -o "$sam_output"
  samtools view -@ 6 -F 0x4 -F 0x2 -bS "$sam_output" > "$OUTPUT_DIR/${assembly_file%.fasta}_$read_file.abyss.bam"
done


# Define spades assembly file names
spades_assembly_files=(
  "ERR204044_trimmed_SPADES_OUT/contigs.fasta"
  "SRR15131330_trimmed_SPADES_OUT/contigs.fasta"
  "SRR18214264_trimmed_SPADES_OUT/contigs.fasta"
)

# Map the original read files to the spades assemblies
for ((i=0; i<${#rspades_assembly_files[@]}; i++)); do
  assembly_file="${spades_assembly_files[i]}"
  read_file="${read_files[i]}"
  sam_output="$OUTPUT_DIR/${read_file}_spades.sam"
  bwa mem -t 6 "$OUTPUT_DIR/$assembly_file" "$INPUT_DIR/${read_file}_1.fastq.gz" "$INPUT_DIR/${read_file}_2.fastq.gz" -o "$sam_output"
  samtools view -@ 6 -F 0x4 -F 0x2 -bS "$sam_output" > "$OUTPUT_DIR/${read_file}_spades.bam"

done

 # Remove .sam files
rm *.sam

# --------------------------
# Next is analysis of mapping results. The mapping fraction aswell as genome coverage from mapped reads.

# Find all _sorted.bam files in the output directory
bam_files=("$OUTPUT_DIR"/*_sorted.bam)

# Loop through the sorted BAM files
for bam_file in "${bam_files[@]}"; do

  # Extract the filename without the extension
  filename=$(basename "$bam_file" _sorted.bam)

  # Calculate the mapping rate using samtools flagstat
  mapping_output=$(samtools flagstat "$bam_file")
  mapped_reads=$(echo "$mapping_output" | awk 'NR==9 {print $1}')
  total_reads=$(echo "$mapping_output" | awk 'NR==1 {print $1}')

  echo "Mapping Output for $filename:"
  echo "$mapping_output"

  # Calculate the mapping rate
  if [[ $total_reads -ne 0 ]]; then
    mapping_rate=$(awk "BEGIN {printf \"%.2f\", ($mapped_reads / $total_reads) * 100}")
  else
    mapping_rate="N/A"
  fi

  echo "Mapping rate for $filename: $mapping_rate%"  >> "$RESULTS_DIR/mapping_results.txt"

  # Calculate the coverage using samtools depth
   if [[ -f "$bam_file" ]]; then
     coverage=$(samtools depth -a "$bam_file" | awk '{ total += $3 } END { printf "%.2f", total / NR }')
     echo "Coverage for $filename: $coverage" >> "$RESULTS_DIR/coverage_results.txt"
   else
     echo "Sorted BAM file not found for $filename" >> "$RESULTS_DIR/coverage_results.txt"
   fi

done

