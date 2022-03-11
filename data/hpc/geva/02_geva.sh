# cat testbatch.txt | awk 'BEGIN{n=1;fn=sprintf("%s.in",n)} NR%2000==0{print fn;n++; fn=sprintf("%s.in",n)}{print > fn}'


# for f in *.in;do
# 	n=${f%.in}
# 	geva_v1beta -t 60 -i test.bin -o run_${n} --positions $f --Ne 100000 --mut 1.2e-8 --hmm ~/.local/geva/hmm/hmm_initial_probs.txt ~/.local/geva/hmm/hmm_emission_probs.txt
# done



while read r;do
	echo $r
	cd chroms/${r}

	geva_v1beta --vcf ${r}.vcf --rec 3.2e-8 --out ${r}_geva

	for f in *.in;do
	 	n=${f%.in}
	 	geva_v1beta -t 60 -i ${r}_geva.bin -o run_${n} --positions $f --Ne 100000 --mut 1.2e-8 --hmm ~/.local/geva/hmm/hmm_initial_probs.txt ~/.local/geva/hmm/hmm_emission_probs.txt
 	done

	cd ../../
done < <(cat input_chrs.txt)
