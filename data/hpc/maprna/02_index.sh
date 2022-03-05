STAR --runThreadN 40 \
	--runMode genomeGenerate \
	--genomeDir ref \
	--genomeFastaFiles GCA_014634065.1_Adig_2.0_genomic.fna \
	--sjdbGTFfile ../../genome/adig-v2-ncbi.gff \
	--sjdbGTFtagExonParentTranscript "Parent" \
	--genomeSAindexNbases 13
