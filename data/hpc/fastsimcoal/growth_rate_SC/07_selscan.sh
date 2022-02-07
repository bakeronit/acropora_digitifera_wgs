cd 3.out.growth_rate_SC

mkdir -p selscan/vcfs
mkdir -p selscan/ihs_out/
mkdir -p selscan/xpehh_out/
mkdir -p selscan/xpnsl_out/

#for f in *.vcf.gz;do
#	echo $f
#	s1=${f#3.out.growth_rate_SC_}
#	s2=${s1%.gen.vcf.gz}
#	for pop in IN NO SO IN_NO IN_SO SO_NO;do
#		for c in $(seq 1 20);do
#			bcftools view -S ${pop}.pop -r contig${c} $f |bgzip > selscan/vcfs/${s2}_${pop}_${c}.vcf.gz
#		done
#	done
#done

## Selscan itself

# iHS scans
# for f in *.vcf.gz;do
# 	echo $f
# 	s1=${f#3.out.growth_rate_SC_}
# 	s2=${s1%.gen.vcf.gz}
# 	for pop in IN NO SO;do
#   		for c in $(seq 1 20);do
# 	  		../../selscan-linux-1.3.0/selscan --ihs --vcf selscan/vcfs/${s2}_${pop}_${c}.vcf.gz --pmap --threads 40 --out selscan/ihs_out/${s2}_${pop}_${c}
#   		done
# 	done
# done


# XP-EHH scans

do_ehh(){
	reflist=$1
 	poplist=$2
  	ref=${reflist%.pop}
  	pop=${poplist%.pop}
	for f in *.vcf.gz;do
		s1=${f#3.out.growth_rate_SC_}
		s2=${s1%.gen.vcf.gz}
  		for c in $(seq 1 20);do
	      	../../selscan-linux-1.3.0/selscan --xpehh --vcf selscan/vcfs/${s2}_${pop}_${c}.vcf.gz \
        	--vcf-ref selscan/vcfs/${s2}_${ref}_${c}.vcf.gz --pmap \
        	--threads 40 --out selscan/xpehh_out/${s2}_${pop}_${c}

			../../selscan-linux-1.3.0/selscan --xpnsl --vcf selscan/vcfs/${s2}_${pop}_${c}.vcf.gz \
        	--vcf-ref selscan/vcfs/${s2}_${ref}_${c}.vcf.gz --pmap \
        	--threads 40 --out selscan/xpnsl_out/${s2}_${pop}_${c}

#        	../../selscan-linux-1.3.0/selscan --ihs --vcf selscan/vcfs/${s2}_${pop}_${c}.vcf.gz --pmap --threads 40 --out selscan/ihs_out/${s2}_${pop}_${c}
	    done
#	    ../../selscan-linux-1.3.0/norm --ihs --files selscan/ihs_out/${s2}_${pop}_*.out --bins 50 --bp-win --min-snps 10 --winsize 50000
	    ../../selscan-linux-1.3.0/norm --ihs --files selscan/xpehh_out/${s2}_${pop}_*.out --bins 50 --bp-win --min-snps 10 --winsize 50000
		../../selscan-linux-1.3.0/norm --ihs --files selscan/xpnsl_out/${s2}_${pop}_*.out --bins 50 --bp-win --min-snps 10 --winsize 50000
  	done
}


do_ehh SO_NO.pop IN.pop
do_ehh IN_NO.pop SO.pop
do_ehh IN_SO.pop NO.pop





# # Norm iHS
# for f in *.vcf.gz;do
# 	echo $f
# 	s1=${f#3.out.growth_rate_SC_}
# 	s2=${s1%.gen.vcf.gz}
# 	for pop in IN NO SO;do
# 	  	../../selscan-linux-1.3.0/norm --ihs --files selscan/ihs_out/${s2}_${pop}_*.out --bins 50 --bp-win --min-snps 10 --winsize 50000
# 	done
# done
