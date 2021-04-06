cat common_uces.tsv | awk '{print $4}' | sort -u > common_probes.txt


# # Extract a list of common probes consisting of exons
rm exon_loci.txt
while read candidate;do
	grep "${candidate}," common_probes.txt >> exon_loci.txt
	grep "${candidate}$" common_probes.txt >> exon_loci.txt
done < <(grep 'design:hexatrans' hexa-v2-sclerac-subset-final-probes.fasta | awk '{print $1}' | sed 's/>//' | sort)

# # And everything else
rm uce_loci.txt
while read candidate;do
	grep "${candidate}," common_probes.txt >> uce_loci.txt
	grep "${candidate}$" common_probes.txt >> uce_loci.txt
done < <(grep '>' hexa-v2-sclerac-subset-final-probes.fasta | grep -v 'design:hexatrans' | awk '{print $1}' | sed 's/>//' | sort)


# # This results in 658 exon loci, 1023 uce loci and 18 shared between them. 
# #

comm -12 <(sort -u exon_loci.txt) <(sort -u uce_loci.txt) > shared_loci.txt


printf "#nexus\nbegin sets;\n" > partition.nex

while read el;do
	on=$(printf "exon_${el%%,*}.fasta")
	pn=$(echo ${el%%,*} | sed 's/-//g' | sed 's/_//g') 
	echo $on
	alignment_len=$(cat alignments/${el}_align.fasta | bioawk -c fastx '{print length($seq)}' | head -n 1)
	cp alignments/${el}_align.fasta iqtree_alignments/${on}

	printf "\tcharset ${pn} = iqtree_alignments/${on}: 1-${alignment_len};\n" >> partition.nex

done < <(comm -13 <(sort -u shared_loci.txt) <(sort -u exon_loci.txt))

while read el;do
	on=$(printf "uce_${el%%,*}.fasta")
	echo $on
	pn=$(echo ${el%%,*} | sed 's/-//g' | sed 's/_//g') 
	alignment_len=$(cat alignments/${el}_align.fasta | bioawk -c fastx '{print length($seq)}' | head -n 1)
	cp alignments/${el}_align.fasta iqtree_alignments/${on}

	printf "\tcharset ${pn} = iqtree_alignments/${on}: 1-${alignment_len};\n" >> partition.nex

done < <(comm -13 <(sort -u shared_loci.txt) <(sort -u uce_loci.txt))

printf "end;\n" >> partition.nex
