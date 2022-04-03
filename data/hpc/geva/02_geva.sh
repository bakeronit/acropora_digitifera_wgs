

while read r;do
	echo $r
	cd chroms/${r}

	geva_v1beta --vcf ${r}.vcf --rec 3.2e-8 --out ${r}_geva

#	rm pos.txt sites.txt
	for f in *.in;do
	 	n=${f%.in}
	 	geva_v1beta -t 60 -i ${r}_geva.bin -o run_${n} --positions $f --Ne 10000 --mut 1.2e-8 --hmm ~/.local/geva/hmm/hmm_initial_probs.txt ~/.local/geva/hmm/hmm_emission_probs.txt
#	 	cat $f >> pos.txt
#	 	cat run_${n} >> sites.txt
 	done

 	


	cd ../../
done < <(cat input_chrs.txt)
