paste variants_snps.tsv <(tail -n+9 output-file-pvalues.txt) | awk -f estsfs2aa.awk > aa.tab

bgzip aa.tab
tabix -s1 -b2 -e2 aa.tab.gz -f

bcftools annotate -a aa.tab.gz -c CHROM,POS,INFO/AA \
	../vcf_file/Adigi.v2.indv74_phased.vcf.gz -h aa.hdr > Adigi.v2.indv74_phased_aa.vcf


cat Adigi.v2.indv74_phased_aa.vcf | awk -f aa2ref.awk > Adigi.v2.indv74_phased_aaref.vcf