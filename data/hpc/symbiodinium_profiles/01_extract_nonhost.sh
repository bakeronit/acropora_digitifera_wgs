
# mkdir flagstats
# for f in $(find /fast/shared/Acropora_digitifera_wgs_bamfile/ -name '*.bam');do
# 	bn=$(basename $f)
# 	sn=${bn%_aligned_duplicates_marked_sorted.bam}
# 	echo $sn
# 	samtools flagstat $f > flagstats/${sn}.flagstat
# done


extract_nh(){
	f=$1
	bn=$(basename $f)
 	sn=${bn%_aligned_duplicates_marked_sorted.bam}

	samtools view -f 13 -b $f > nonhostbam/${sn}_nonhost.bam
}

export -f extract_nh

parallel -j 40 extract_nh ::: $(find /fast/shared/Acropora_digitifera_wgs_bamfile/ -name '*.bam' | tr '\n' ' ')

cd nonhostbam

for f in *_nonhost.bam;do
	echo $f
	bn=${f%_nonhost.bam}
	samtools fastq -1 ${bn}_R1.fastq -2 ${bn}_R2.fastq -s ${bn}_s.fastq $f
done

mkdir ../reads
mv *.fastq ../reads/

