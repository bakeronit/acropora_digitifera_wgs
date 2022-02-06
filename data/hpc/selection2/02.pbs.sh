mkdir -p pbs

f=../vcf_file/Adigi.v2.filtered.vcf.gz

plink2 --vcf $f --fst site report-variants --pheno populations.txt --allow-extra-chr --out pbs/plink2

echo "IN NO SO" > pbs/plink2.pbs
paste pbs/plink2.IN.NO.fst.var pbs/plink2.IN.SO.fst.var pbs/plink2.NO.SO.fst.var | awk -f plinkfst2pbs.awk >> pbs/plink2.pbs


cat pbs/plink2_noheader.pbs | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4,$5}' > pbs/plink2_noheader.bed
bedtools intersect -b ../tracks/sweeps.gff3 -a pbs/plink2_noheader.bed > pbs/pbs_sweeps.tsv

bedtools subtract -a pbs/plink2_noheader.bed -b pbs/pbs_sweeps.tsv > pbs/pbs_nosweeps.tsv

bedtools sort -i ../tracks/sweeps.gff3 > sweeps_sorted.gff
bedtools sort -i pbs/plink2_noheader.bed > pbs/plink2_noheader_sorted.bed

bedtools map -a sweeps_sorted.gff -b pbs/plink2_noheader_sorted.bed -c 4 -o max > pbs/max_pbs_in_sweeps.tsv