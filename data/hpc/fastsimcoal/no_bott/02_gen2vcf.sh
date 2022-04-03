cd 3.out.no_bott
for f in *.gen;do
	echo $f
	if [[ ! -f ${f}.vcf.gz ]]; then
		cat ${f} | awk -f ../../gen2vcf.awk > ${f}.vcf
		bgzip ${f}.vcf
	fi
done

