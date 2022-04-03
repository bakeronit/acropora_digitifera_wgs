# We just want broad summaries of fst here to compare with the bottleneck scenario

cd 3.out.no_bott

mkdir pca

for f in *.vcf.gz;do
	echo $f
	s1=${f#3.out.growth_rate_SC_}
	s2=${s1%.gen.vcf.gz}

	plink2 --vcf ${f} --pca --pheno phenotypes.txt --allow-extra-chr --out pca/${s2}

done
