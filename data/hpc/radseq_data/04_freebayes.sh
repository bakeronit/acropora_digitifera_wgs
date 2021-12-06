# Set strict input base call and mapping qualities

for dataset in allrad;do
	freebayes-parallel genome_regions.txt 46 -f GCA_014634065.1_Adig_2.0_genomic.fna \
		-L ${dataset}_bamlist.txt \
		--genotype-qualities \
		-E -1 \
		-m 30 -q 20 \
		-K --strict-vcf > ${dataset}.vcf
done


