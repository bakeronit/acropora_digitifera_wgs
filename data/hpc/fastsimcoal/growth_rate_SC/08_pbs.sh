cd 3.out.growth_rate_SC

mkdir -p pbs

#for f in *.vcf.gz;do
for f in 3.out.growth_rate_SC_16_1.gen.vcf.gz;do	
	echo $f
	s1=${f#3.out.growth_rate_SC_}
	s2=${s1%.gen.vcf.gz}

	plink2 --vcf $f --fst site report-variants --pheno phenotypes.txt --allow-extra-chr --out pbs/${s2}

	echo "IN NO SO" > pbs/${s2}.pbs
	paste pbs/${s2}.IN.NO.fst.var pbs/${s2}.IN.SO.fst.var pbs/${s2}.NO.SO.fst.var | awk -f ../../plinkfst2pbs.awk >> pbs/${s2}.pbs

done

