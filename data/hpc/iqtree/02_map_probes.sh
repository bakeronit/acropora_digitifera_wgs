# For each of our target genomes we need to map the Hexa probeset
#
# Use samclip to retain alignments only if they are clipped by less than 20 bp
#
module load bwa/0.7.17

bwa index aten_final_0.11.fasta
bwa index GCA_014634065.1_Adig_2.0_genomic.fna
bwa index GCF_004143615.1_amil_sf_1.1_genomic.fna


bwa mem -M aten_final_0.11.fasta hexa-v2-sclerac-subset-final-probes.fasta > aten_probes.sam
samclip --ref aten_final_0.11.fasta --max 20 < aten_probes.sam > aten_probes_clip.sam


bwa mem -M GCA_014634065.1_Adig_2.0_genomic.fna hexa-v2-sclerac-subset-final-probes.fasta > adi_probes.sam
samclip --ref GCA_014634065.1_Adig_2.0_genomic.fna --max 20 < adi_probes.sam > adi_probes_clip.sam

bwa mem -M GCF_004143615.1_amil_sf_1.1_genomic.fna hexa-v2-sclerac-subset-final-probes.fasta > amil_probes.sam
samclip --ref GCF_004143615.1_amil_sf_1.1_genomic.fna --max 20 < amil_probes.sam > amil_probes_clip.sam
