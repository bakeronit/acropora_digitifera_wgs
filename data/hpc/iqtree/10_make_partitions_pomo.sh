printf "#nexus\nbegin sets;\n" > pomo_partition.nex

while read el;do
	on=$(printf "exon_${el%%,*}.fasta")
	pn=$(echo ${el%%,*} | sed 's/-//g' | sed 's/_//g') 
	echo $on
	alignment_len=$(cat alignments/${el}_align.fasta | bioawk -c fastx '{print length($seq)}' | head -n 1)
	cat pomo_select.txt | xargs samtools faidx ../pomo/alignments/${el}_align.fasta > pomo_alignments/${on}
#	cp alignments/${el}_align.fasta pomo_alignments/${on}

	printf "\tcharset ${pn} = pomo_alignments/${on}: 1-${alignment_len};\n" >> pomo_partition.nex

done < <(comm -13 <(sort -u shared_loci.txt) <(sort -u exon_loci.txt))

while read el;do
	on=$(printf "uce_${el%%,*}.fasta")
	echo $on
	pn=$(echo ${el%%,*} | sed 's/-//g' | sed 's/_//g') 
	alignment_len=$(cat alignments/${el}_align.fasta | bioawk -c fastx '{print length($seq)}' | head -n 1)

	cat pomo_select.txt | xargs samtools faidx ../pomo/alignments/${el}_align.fasta > pomo_alignments/${on}
#	cp alignments/${el}_align.fasta pomo_alignments/${on}

	printf "\tcharset ${pn} = pomo_alignments/${on}: 1-${alignment_len};\n" >> pomo_partition.nex

done < <(comm -13 <(sort -u shared_loci.txt) <(sort -u uce_loci.txt))

printf "end;\n" >> pomo_partition.nex
