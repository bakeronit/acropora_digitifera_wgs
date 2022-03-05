for f in *_sorted.bam;do
	echo $f;
	samtools view $f 'BLFC01000154.1' -b > ${f%_sorted.bam}_BLFC01000154.1.bam
done

samtools merge BLFC01000154.1.bam *_BLF*.bam -O bam -f
samtools index BLFC01000154.1.bam