for f in *out.bam;do 
	sample=${f%Aligned*};
	echo $sample;
	samtools sort -@ 10 -O bam -o ${sample}_sorted.bam $f 
	samtools index ${sample}_sorted.bam
done
