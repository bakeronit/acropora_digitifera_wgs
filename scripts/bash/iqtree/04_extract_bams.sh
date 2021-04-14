module load bedtools2/2.25.0

# Extract data from bam files

bedtools sort -faidx aten_final_0.11.fasta.fai -i aten_probes_clip.bed > aten_probes_sorted.bed
bedtools merge -i aten_probes_sorted.bed -c 4 -o collapse > aten_probes_merged.bed

bedtools sort -faidx GCA_014634065.1_Adig_2.0_genomic.fna.fai -i adi_probes_clip.bed > adi_probes_sorted.bed
bedtools merge -i adi_probes_sorted.bed -c 4 -o collapse > adi_probes_merged.bed

bedtools sort -faidx GCF_004143615.1_amil_sf_1.1_genomic.fna.fai -i amil_probes_clip.bed > amil_probes_sorted.bed
bedtools merge -i amil_probes_sorted.bed -c 4 -o collapse > amil_probes_merged.bed







# for f in *merged_marked.bam;do
# 	echo $f
# 	samtools index $f
# 	samtools view -M -b -L aten_probes_merged.bed $f > ${f%.bam}_regions.bam
# done


while read f;do
	ln -s $f .
	fb=$(basename $f)
	echo $fb
	samtools index $fb
	samtools view -M -b -L adi_probes_merged.bed $fb > ${fb%.bam}_regions.bam
done < adi_bams.txt

