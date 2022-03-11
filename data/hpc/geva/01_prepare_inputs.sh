# Estimating SNP ages using GEVA
# https://github.com/pkalbers/geva

# Export chromosomes individually and convert to numeric ids
#
mkdir -p chroms

while read r;do
	echo $r
	mkdir -p chroms/${r}/
	bcftools view ../vcf_file/Adigi.v2.indv74_phased.vcf.gz --regions ${r} | \
		sed s/${r}/1/ >  chroms/${r}/${r}.vcf

	cd chroms/${r}

	cat ${r}.vcf | awk '!/^#/{print $2}' | awk 'BEGIN{n=1;fn=sprintf("%s.in",n)} NR%2000==0{print fn;n++; fn=sprintf("%s.in",n)}{print > fn}'
	cd ../../
done < <(cat input_chrs.txt)
