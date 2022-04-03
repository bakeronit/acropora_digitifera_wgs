# Info to print out
# ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
# ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">

# Note. POS0 is used so these positions are in 0-based coordinates (Like BED)

#bcftools query ../vcf_file/Adigi.v2.filtered.vcf.gz -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT{0}\t%AC\t%AN\n" > variants.bed


bedtools intersect \
	-a variants.bed \
	-b <(cat snps.tsv | grep -v '^ref' | sort -k1,1 -k2,2n  | awk 'NF==5{printf("%s\t%s\t%s\t%s,%s,%s\n",$1,$2,$2+1,$3,$4,$5,$6)}') \
	-loj -F 1 -f 1 > variants_snps.tsv



cat variants_snps.tsv | awk -f snps2estsfs.awk > estsfs-input.txt