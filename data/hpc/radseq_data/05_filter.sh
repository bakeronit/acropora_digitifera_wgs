# Filter for sites with high MAF and high quality genotype calls. 
# These sites should be ideally suited to relatedness statistics

vcftools --vcf allrad.vcf --maf 0.1 \
	--min-meanDP 15 --hwe 0.001 \
	--minDP 8 --min-alleles 2 --max-alleles 2 \
	--recode --stdout > allrad_filt.vcf


bcftools merge  Adigi.v2.filtered.vcf.gz allrad_filt.vcf.gz | 
	bcftools view -i 'F_MISSING<0.1' > allrad_merged.vcf 