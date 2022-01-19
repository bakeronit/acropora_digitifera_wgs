#cat allgo_genes.tsv | awk -F '\t' '{printf("%s.t1\n", $5)}' | grep -v 'genes' | xargs samtools faidx protein.fa > allgo_genes_t1.fasta

#cat allgo_genes_t1.fasta | bioawk -c fastx '{print $name,$seq}' | sed 's/.t1//' > allgo_genes_t1.tsv


#~/bin/ncbi-blast-2.12.0+/bin/blastp -db nr -query allgo_genes_t1.fasta -remote -entrez_query 'NOT Acropora digitifera [ORGN]' -outfmt "6 std stitle ssciname staxid " -max_hsps 1 -evalue 0.0001 >  allgo_genes_t1.blastp

cat allgo_genes_t1.fasta | bioawk -c fastx 'BEGIN{b=0;fname=sprintf("batch_%s.fasta",b)} NR%3==0{b+=1;fname=sprintf("batch_%s.fasta",b)} {printf(">%s\n%s\n",$name,$seq) >> fname}'

for qf in batch*.fasta;do

	echo "Blasting batch ${qf}"

	~/bin/ncbi-blast-2.12.0+/bin/blastp -db nr -query $qf -remote -entrez_query 'NOT Acropora digitifera [ORGN]' -outfmt "6 std stitle ssciname staxid " -max_hsps 1 -evalue 0.0001 >  ${qf%.fasta}.blastp

	sleep 100

done