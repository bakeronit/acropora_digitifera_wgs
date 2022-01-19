cd 3.out.growth_rate_SC

mkdir -p admixture

cd admixture

for f in ../*.vcf.gz;do
    echo $f
    sf=${f#../3.out.growth_rate_SC_}
    sample=${sf%.gen.vcf.gz}

	gunzip -c $f | sed s/contig// > ${sample}.vcf
	~/bin/plink2 --vcf ${sample}.vcf -out ${sample} --make-bed --indep-pairwise 50 10 0.1 --set-missing-var-ids @:#
	~/bin/plink2 --bfile ${sample} --extract ${sample}.prune.in --allow-extra-chr --out ${sample}_prune --make-bed
	../../../admixture_linux-1.3.0/admixture ${sample}_prune.bed 3 -j10
done
