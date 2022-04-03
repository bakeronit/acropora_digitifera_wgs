mkdir -p regions

while read line;do
	parts=( $line )
	sweep_id=${parts[0]}
	position=${parts[1]}
	allele=${parts[2]}
	chrom=${parts[3]}
	start=${parts[4]}
	end=${parts[5]}

	echo $chrom $position $allele ${sweep_id}

	vcf_path="regions/${sweep_id}_aaref.vcf"

	# Extract vcf for the region
	if [[ ! -f ${vcf_path} ]]; then
		echo "Exporting vcf file ${chrom} $start $end"
		vcftools --gzvcf ../geva/Adigi.v2.indv74_phased_aaref.vcf.gz --recode --recode-INFO-all --stdout --chr $chrom --from-bp $start --to-bp $end  > ${vcf_path}
	fi

done < <(cat sweep_sites.tsv | grep -v 'sweep_id')

