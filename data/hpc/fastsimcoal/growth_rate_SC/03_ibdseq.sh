cd 3.out.growth_rate_SC

mkdir -p ibdseq

#for f in *.vcf.gz;do
for in 3.out.growth_rate_SC_16_1.gen.vcf.gz;do
	echo $f
	for c in $(seq 1 100);do
		s1=${f#3.out.growth_rate_SC_}
		s2=${s1%.gen.vcf.gz}
		java -Xmx2000m -jar ../../../ibdseq/ibdseq.r1206.jar gt=${f} out=ibdseq/${s2}.${c} nthreads=4 chrom=contig${c}
	done
done