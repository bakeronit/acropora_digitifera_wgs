for f in *_1.fastq;do 
	sample=${f%_1.fastq}
	echo $sample
	STAR --runThreadN 60 \
	--genomeDir ref \
	--readFilesIn ${sample}_1.fastq ${sample}_2.fastq \
	--outFileNamePrefix ${sample} \
	--outSAMtype BAM Unsorted \
	--outSAMunmapped Within \
	--outSAMattributes Standard \
	--quantMode TranscriptomeSAM GeneCounts
done
