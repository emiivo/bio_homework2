#!/bin/bash

# Download raw sequencing reads
prefetch -O ~/HW2/inputs ERR204044 SRR15131330 SRR18214264

# Define input directory and files
INPUT_DIR=~/HW2/inputs
OUTPUT_DIR=~/HW2/outputs
FILES=(ERR204044 SRR15131330 SRR18214264)

# Loop through input files and run fastq-dump and fastqc for each file
for file in "${FILES[@]}"
do
  fastq-dump --outdir "$INPUT_DIR" --gzip --split-files "$INPUT_DIR/$file"
  fastqc "$INPUT_DIR/${file}_1.fastq.gz" "$INPUT_DIR/${file}_2.fastq.gz" -o ~/HW2/outputs/
done

# The quality of the sequences is ok, but not great, there are some adapters, there is some sequence duplication.

for file in "${FILES[@]}"
do
  # Define input and output file names
  INPUT1="${INPUT_DIR}/${file}_1.fastq.gz"
  INPUT2="${INPUT_DIR}/${file}_2.fastq.gz"
  BASENAME="${file}_trimmed"


  # Run trim_galore
  trim_galore --paired --illumina --phred33 --quality 20 --length 20 --paired ${INPUT1} ${INPUT2} --output_dir ${OUTPUT_DIR} --basename ${BASENAME}

done

# Run FastQC for trimmed files
for file in "${FILES[@]}"
do
  fastqc "$OUTPUT_DIR/${file}_trimmed_val_1.fq.gz" "$OUTPUT_DIR/${file}_trimmed_val_2.fq.gz" -o ~/HW2/outputs/
  mv "$OUTPUT_DIR/${file}_trimmed_val_1.fq.gz" "$OUTPUT_DIR/${file}_trimmed_1.fastq.gz"
  mv "$OUTPUT_DIR/${file}_trimmed_val_2.fq.gz" "$OUTPUT_DIR/${file}_trimmed_2.fastq.gz"
  fastqc "$OUTPUT_DIR/${file}_trimmed_1.fastq.gz" "$OUTPUT_DIR/${file}_trimmed_2.fastq.gz" -o ~/HW2/outputs/

done

# The quality of the files has improved, around 1% of basepairs have been filtered out, overall data looks better. Sequence duplication levels are sill
# not great, but everything else looks ok. Proportion per base sequence content is not very good.

multiqc ~/HW2/outputs/ -o ~/HW2/outputs/

