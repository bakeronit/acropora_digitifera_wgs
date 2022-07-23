# Our goal is to place A. digitifera in WA and Japan within a broader evolutionary context.
# For this we will build a tree based on UCE sequences present in;
# A. millepora
# A. digitifera - Japan
# A. digitifera - North offshore, inshore, South offshore
# A. tenuis - Magnetic Island
# A. tenuis - Fitzroy Island

# We will need to
# 0. Reference genomes for all three species
# 1. Identify the genomic locations of UCE's in all three species reference genomes
# 2. Call consensus sequences from bams for WGS data
# 3. Extract reference sequences in other cases

# Probe-set
# This is the v2 probeset published by Cowman et al 10.1016/j.ympev.2020.106944


# Reference Genomes
#
# Amil v1 (NCBI)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/143/615/GCF_004143615.1_amil_sf_1.1/GCF_004143615.1_amil_sf_1.1_genomic.fna.gz
# Adig v2 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/634/065/GCA_014634065.1_Adig_2.0/GCA_014634065.1_Adig_2.0_genomic.fna.gz
# Aten
wget http://aten.reefgenomics.org/download/aten_final_0.11.fasta.gz


# Bam files
#
# Aten MI
cp ~/Projects/atenuis_wgs/hpc/gatk3/MI-1-4_merged_marked.bam .
# Aten FI
cp ~/Projects/atenuis_wgs/hpc/gatk3/FI-1-3_merged_marked.bam .

# Adi JP
#
#../pcangsd/jp_bam.txt 
#
# Adi WA
#../pcangsd/wa_bam.txt  
#



