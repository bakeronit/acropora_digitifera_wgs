# bwa 0.7.17-r1188
# samtools 1.14
# 
map_reads(){

	mitogenome=$1

	bamfile=$2

	if [ ! -f ${sample}_mitoreads.bam ];then

		sample=${bamfile%_aligned_duplicates_marked_sorted.bam}

		samtools fastq -F 1024 ${bamfile} | bwa mem -t 16 ${mitogenome} - | samtools view -b -F 4 - > ${sample}_mitoreads.bam
	fi

}

export -f map_reads


parallel -j 40 map_reads 'NC_022830.1.fasta' ::: $(ls *_aligned_duplicates_marked_sorted.bam | tr '\n' ' ')