mkdir -p regions

while read line;do

	parts=( $line )
	sweep_id=${parts[0]}
	position=${parts[1]}
	allele=${parts[2]}
	chrom=${parts[3]}
	start=${parts[4]}
	end=${parts[5]}
	pop=${parts[6]}


	snp_id=${sweep_id}_${position}

	echo $position $allele ${sweep_id}

	vcf_path="regions/${sweep_id}_aaref.vcf"

	# Extract vcf for the region
	if [[ ! -f ${vcf_path} ]]; then
		echo "No vcf at path ${vcf_path}"
		exit
	fi

	for rep in $(seq 1 10);do
		echo "Rscript run_startmrca.R $position $allele ${vcf_path} ${snp_id} $pop $rep"

		cat qsub_template.sh | sed "s/startmrca/${rep}_${snp_id}/" > ${snp_id}_${rep}.qsub
		echo "Rscript run_startmrca.R $position $allele ${vcf_path} ${snp_id} $pop $rep" >> ${snp_id}_${rep}.qsub
		qsub ${snp_id}_${rep}.qsub
	done

done < <(cat sweep_sites.tsv | grep -v 'sweep_id')

