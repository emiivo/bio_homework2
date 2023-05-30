#!/bin/bash

REFERENCE_DIR=~/HW2/references
OUTPUT_DIR=~/HW2/outputs

# Download reference files for CP015498 nucleotides and proteins:

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&dopt=fasta&sendto=on&id=CP015498.1" -O "$REFERENCE_DIR/CP015498.fasta"

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=1318702381&db=nuccore&report=fasta&retmode=text" -O "$REFERENCE_DIR/CP015498_protein.fasta"

# ---------------------------------------
# First step for genome analysis and annotation is using Gepard to create dotplots between samples. I did this through my own computer, 
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

# --------------------------------------
# Next I analysed my genomes with BUSCO. I also did this through my own computer, through galaxy.eu.
# The results of this analysis is further proving my point that there is something wrong with the ABySS assemblies, I uploaded the scaffolded data of all
# six datasets (3 from ABySS and 3 from SPAdes). 

# Firstly, the ABySS assemblies:
# 1) ERR204044 - overall, the assembly is pretty good, 98.6% of the BUSCO groups were identified.
# Out of the 402 BUSCO groups searched, 396 were complete, with 360 being single-copy and 36 being duplicated. 
# Only 1% of the BUSCOs were fragmented, indicating that most of the identified genes were complete.

# 2) SRR15131330 - many of the BUSCOs (15.2%) were fragmented, and most were not even found in the assembly (80.1%).
# Only 4.7% of the BUSCOs in the lineage dataset were identified in the assembly. Overall, assembly failed the analysis.

# 3) SRR18214264 - again, the majority of BUSCOs (60.4%) were not found in the assembly.
# 21.4% of the BUSCOs were identified as fragmented, indicating partial gene sequences and only 18.2% of the BUSCOs were identified.

# Overall, results of ABySS assemblies are not good.

# Next, the SPAdes assemblies:
# 1) ERR204044 - 99.0% of the BUSCOs in the lineage dataset were identified in the assembly. Only 0.5% of the BUSCOs were identified as fragmented, 
# indicating nearly complete gene sequences and 0.5% were missing. This result is very good.
# 2) SRR15131330 - 98.5% of the BUSCOs were identified in the assembly. Only 0.5% of the BUSCOs were identified as fragmented and 
# only 1 of the BUSCOs were missing from the assembly. This result is also good.
# 3) SRR18214264 - ad with ERR204044, 99.0% of the BUSCOs in the lineage dataset were identified in the assembly. Only 0.5% of the BUSCOs were identified as fragmented, 
# indicating nearly complete gene sequences and 0.5% were missing. This result also very good.

# --------------------------------------

# Next step is gene prediction with GeneMarkS-2. This was done on my own computer through the graphical user interface.
# Results will be discussed when I compare all three of the required gene prediction tools.

# --------------------------------------

# The next step is to using RAST genome annotation server, predict and annotate genes in my assemblies. 
# I used RAST in my own computer as well. The results will be discussed in a later part of the code, but a problem I encountered with my assemblies 
# made with ABySS are that SRR18214264 and SRR15131330 were too low quality to use 
# (Message of the site - "We are sorry, but your job appears to still be too low quality for RAST to analyze").
# I decided to still do the rest of the steps with them because I hope that maybe I will still get at least some useful result or get some idea of why
# exactly their assembly failed.

# --------------------------------------
# Next step - gene prediction using CP015498 genes and proteins models as well as local blast.

CP015498_NUCLEOTIDE=$REFERENCE_DIR/CP015498.fasta
CP015498_PROTEIN=$REFERENCE_DIR/CP015498_protein.fasta

# Set the path to the assembly files (scaffolded)
ASSEMBLY_DIRECTORIES=(
    "$OUTPUT_DIR/ERR204044_spades_reordered"
    "$OUTPUT_DIR/ERR204044_abyss_reordered"
    "$OUTPUT_DIR/SRR15131330_spades_reordered"
    "$OUTPUT_DIR/SRR15131330_abyss_reordered"
    "$OUTPUT_DIR/SRR18214264_abyss_reordered"
    "$OUTPUT_DIR/SRR18214264_spades_reordered"
)

# Set the path for the BLAST analysis directory
BLAST_DIR="$OUTPUT_DIR/blast_analysis"

# Create the BLAST analysis directory if it doesn't exist
mkdir -p "$BLAST_DIR"

# Create BLAST database from the reference nucleotide sequence file
makeblastdb -in "$CP015498_NUCLEOTIDE" -dbtype nucl -out "$BLAST_DIR/nucleotide"

# Create BLAST database from the reference protein sequence file
makeblastdb -in "$CP015498_PROTEIN" -dbtype prot -out "$BLAST_DIR/protein"

# Perform gene prediction for each assembly directory
for ASSEMBLY_DIR in "${ASSEMBLY_DIRECTORIES[@]}"; do
    # Extract the assembly directory name without the path and remove the word "reordered"
    ASSEMBLY_NAME=$(basename "$ASSEMBLY_DIR" | sed 's/_reordered$//')

    # Get the assembly file path
    ASSEMBLY_FILE="$ASSEMBLY_DIR/ragtag.scaffold.fasta"

    # Run blastn
    blastn -query "$ASSEMBLY_FILE" -db "$BLAST_DIR/nucleotide" -evalue 1e-6 -num_threads 4 -out "$BLAST_DIR/${ASSEMBLY_NAME}_blastn_results.txt"

    # Run blastx
    blastx -query "$ASSEMBLY_FILE" -db "$BLAST_DIR/protein" -evalue 1e-6 -num_threads 4 -out "$BLAST_DIR/${ASSEMBLY_NAME}_blastx_results.txt"
    
done

# ------------------------------------------
# Now I have all three required types of gene predictions - made with GeneMarkS-2, RAST and local BLAST. I will be comparing them.

# For GeneMarkS-2, I calculated genes and overlaps manually using LibreOffice Calc. I got overlapping sequences by 
# comparing where the prior sequence ends and the next one begins, if the sequence begins in a position where the sequence before
# hasn't ended, I considered them overlapping.

# For RAST I used the same method. Note that there are no gene predictions for assemblies SRR15131330 and SRR18214264 made with ABySS.

# 1) ERR204044 assembly made with SPAdes - 1945 hits, 219 overlap; (GENEMARKS), 2683 hits, 316 overlap (RAST)

# 2) ERR204044 assembly made with ABySS - 1945 hits, 218 overlap (GENEMARKS); 3454 hits, 622 overlap (RAST)

# 3) SRR15131330 assembly made with SPAdes - 2024 hits, 210 overlaps (GENEMARKS); 2854 hits, 336 overlap (RAST)

# 4) SRR15131330 assembly made with ABySS - 3666 hits - all marked as atypical, 36 overlaps (GENEMARKS). I believe this data is not trustworthy.

# 5) SRR18214264 assembly made with SPAdes - 1944 hits, 218 overlaps (GENEMARKS); 2622 hits, 268 overlap (RAST)

# 6) SRR18214264 assembly made with ABySS -  3365 hits, 107 overlaps all ATPYICAL (GENEMARKS). I believe this data is not trustworthy.

# For the next steps I decided to only continue analysing the SPAdes assembly.

# ----------------------------------------

# Next step - create a ring diagram using RAST predictions. I did this on the RAST page, in seed viewer, by choosing sequence based comparisson 
# and comparing all three predictions to reference Lactobacillus helveticus DPC 4571 (already in RAST database).
# The results show a high similarity to the reference for all predicions (mostly from 95% to 100%). The most similar sequences to each other are 
# ERR204044 and SRR18214264, and the most similar sequence to the reference is SRR15131330.

# ----------------------------------------

# Next I created a phylogenetic tree using 16S sequences including the reference genome CP015498 and an outgroup (AP017922 - Staphylococcus aureus).
# I got all the 16S genes by downloading gff3 format and FASTA of the reference and outgroup, downloading the coordinates of these genes in my assemblies 
# by filtering the annotations created by RAST in the RAST Seed Viewer. Next using LibreOffice Calc I filtered out the 16S genes from reference and outgroup
# and used bedtools (in my own computer - I decided it would be simpler because I have all the files already on it) to extract the sequences of the specific coordinates
# from FASTA files (for my ssemblies - scaffolded). As an example, I used a similar command for all my files:

# bedtools getfasta -fi CP015498.fasta -bed CP015498_16S.bed -fo reference_16S_sequences.fasta

# Then I used Jalview to upload all the sequences, I alligned the sequences with the Muscle default settings because Muscle is an aligning tool that is
# highly accurate. After aligning, I created the tree using the options for DNA and Neighbour Joining. I chose the Neighbor Joining method instead of 
# Average Distance because Neighbor joining is efficiencient and handles large datasets well.

# -----------------------------------------

# Next I selected five proteins and created a multi-gene tree. Again, I used the same reference and outgroup and used my BUSCO predictions for this step.
# My proteins - DNA polymerase III PolC-type, DNA-directed RNA polymerase subunit beta, DNA-directed RNA polymerase subunit beta',
# helicase-exonuclease AddAB subunit AddA, valine--tRNA ligase. 

# Again, I used the same method to format them so I can use them in Jalview. 
# This time to create the tree I chose options Neighbour Joining and PID. PID (Percent Identity) is a useful metric for comparing proteins because it 
# quantifies the degree of sequence similarity, indicating the level of conservation and potential functional similarity between protein sequences.

# -----------------------------------------

# My phylogenetic trees look similar and share clusters. This is expected because I am using the same genomes for comparison.
# The presence of shared clusters in the phylogenetic trees supports the expected evolutionary relationships among the genomes,
# indicating common ancestry and similar evolutionary trajectories.
# The clusters that stood out in the tree, as expected, were of the outgroup as it is a different organism than the rest.

# -----------------------------------------

# Looking at my Gepard comparisons, ring diagram and phylogenetic trees, I believe that the assemblies of ERR204044 and SRR18214264 are the most similar to
# each other. In the Gepard dotplot, even though the plot is not completely even, it shows the most even negative coorelation, which means these sequences
# align well with each other and share the most similarities. This is also seen in the ring diagram, because when comparing them to the reference, they 
# share very similar results. The same trend exists in the phylogenetic trees, as these two assemblies can be seen mostly clustered together.


