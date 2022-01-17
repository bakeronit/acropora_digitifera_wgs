for f in $(find /fast/shared/Acropora_digitifera_wgs_bamfile/ -name '*.bam');do
	echo $f
	bn=$(basename $f)
	sn=${bn%_aligned_duplicates_marked_sorted.bam}
	samtools fastq -1 reads_all/${sn}_R1.fastq.gz -2 reads_all/${sn}_R2.fastq.gz $f
done

