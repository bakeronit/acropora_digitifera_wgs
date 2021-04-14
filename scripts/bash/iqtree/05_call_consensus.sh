# Extract data from bam files

module load bcftools/1.11
module load parallel/20180122

do_consensus(){
	ref=$1
	aln=$2
	calls=${aln%.bam}

	bcftools mpileup -Ou -f $ref $aln | \
	bcftools call -mv -Oz -o ${calls}.vcf.gz

	bcftools index ${calls}.vcf.gz


	# normalize indels
	bcftools norm -f ${ref} ${calls}.vcf.gz -Ob -o ${calls}.norm.bcf

	# filter adjacent indels within 5bp
	bcftools filter --IndelGap 5 ${calls}.norm.bcf -Ob -o ${calls}.norm.flt-indels.bcf

	bcftools index ${calls}.norm.flt-indels.bcf

	bcftools consensus -I -f ${ref} ${calls}.norm.flt-indels.bcf > ${calls}.consensus.fa

}

export -f do_consensus

# A tenuis
parallel -J 2 do_consensus aten_final_0.11.fasta ::: FI-1-3_merged_marked_regions.bam MI-1-4_merged_marked_regions.bam

# A digitifera
parallel -J 15 do_consensus GCA_014634065.1_Adig_2.0_genomic.fna ::: *aligned*regions.bam



