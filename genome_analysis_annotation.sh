
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
