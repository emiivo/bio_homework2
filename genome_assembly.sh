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

# ---------------------------------
# Mapping rate: 
# For the ERR204044 dataset, the mapping rate for the ABySS assembly is 86.28%, while for the Spades assembly, it is 98.83%.
# The SPAdes assembly has a higher mapping rate, it aligns more reads to the reference genome compared to the ABySS assembly.
# This suggests that the SPAdes assembly produced a more complete and accurate representation of the genome.

# For the SRR15131330 dataset, the mapping rate for the ABySS assembly is 48.40%, while for the SPAdes assembly, it is 98.72%. 
# Again, the SPAdes assembly exhibits a significantly higher mapping rate, implying a better alignment of reads and potentially 
# a more reliable assembly compared to the ABySS assembly.

# Lastly, for the SRR18214264 dataset, the mapping rate for the ABySS assembly is 46.93%, while for the SPAdes assembly, it is 98.18%.
# The Spades assembly demonstrates a higher mapping rate, indicating a superior alignment of reads.

# In summary, the mapping rates for the SPAdes assemblies consistently outperform those of the ABySS assemblies in all three datasets provided.
# This is not exactly the result I was expecting, because when using QUAST, all assemblies seemed to be similar in quality, in ERR204044 
# it seemed that ABySS generated a better assembly.
# ---------------------------------
# Genome coverage:
# For the ERR204044 dataset, the coverage for the Abyss assembly is 29.80, while for the Spades assembly, it is 19.95. 
# The Abyss assembly has a higher coverage, suggesting that it may have captured more sequencing depth 
# and potentially provided a more comprehensive representation of the genome compared to the Spades assembly.

# For the SRR15131330 dataset, the coverage for the Abyss assembly is 434.17, while for the Spades assembly,
# it is 218.10. The Abyss assembly demonstrates a higher coverage value, indicating a greater sequencing depth.

# Lastly, for the SRR18214264 dataset, the coverage for the Abyss assembly is 42.80, while for the Spades assembly,
# it is 6.56. Once again, the Abyss assembly exhibits a higher coverage value.

# The mapping rates are higher for the SPAdes assemblies, indicating that a larger percentage 
# of reads were successfully aligned to the reference genome using the SPAdes method. However, the coverage values are
# higher for ABySSs assemblies, indicating a higher sequencing depth.


# ---------------------------------------
# Next step is using Gepard to create dotplots between samples. I did this through my own computer, 
# results can be found in moodle.
# I compared both SPAdes and ABySS scaffolds, but the I analyse AbySS assemblies, the more I believe my ABySS assemblies have either been made incorrectly
# or I compromised them in a different step. The low mapping could an indicator of the assemblies not being correct. The dotplots do show coorelation,
# but they are distorted and I am not sure they can even be analised properly.

# As for the SPAdes samples, these are the results I got:
# 1) ERR204044 and SRR15131330 - the dotplot indicates the sequences align in many areas, suggesting the sequences are similar.
# This suggests that they may share common genomic regions or have similar genetic content.
# 2) SRR18214264 and ERR204044 - the dotplot also shows that the sequences are similar, but they have more areas that differ in comparrisson to 
# the ERR204044 and SRR15131330 dotplot. The differences observed in the dotplot suggest variations or divergent regions between the two samples,
# indicating distinct genomic features or genetic differences.
# 3) SRR18214264 and SRR15131330 - the dotplot is very similar to dotplot of ERR204044 and SRR15131330, indicating a high level of sequence 
# alignment and similarity between these two samples. 
