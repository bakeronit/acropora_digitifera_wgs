#conda activate
# Bedtools v2.30.0

for pop in inshore north_offshore south_offshore;do
	bl="${pop}_bam.txt"
#	bedtools genomecov -bg -ibam -i $(cat ${bl} | awk '{printf("-i %s ",$1)}') -g GCA_014634065.1_Adig_2.0_genomic.fna > ${pop}_readcov.bg

	cat ${pop}_readcov.bg | sort -k1,1 -k2,2n > ${pop}_readcov_sorted.bg 
	bedGraphToBigWig ${pop}_readcov_sorted.bg adig2_chromsizes.txt ${pop}_readcov.bw
done