cd 3.out.growth_rate_SC

mkdir -p plink2

#for f in *.vcf.gz;do
for f in 3.out.growth_rate_SC_16_1.gen.vcf.gz;do
	echo $f
	s1=${f#3.out.growth_rate_SC_}
	s2=${s1%.gen.vcf.gz}
	~/bin/plink2 --vcf ${f} --allow-extra-chr --het
	mv plink2.het plink2/${s2}.het
done

