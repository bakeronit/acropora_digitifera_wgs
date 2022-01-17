
call_cons(){

	mitogenome=$1

	bamfile=$2

	sample=${bamfile%_mitoreads.bam}

	samtools sort ${bamfile} | \
	samtools mpileup -uf ${mitogenome} - | \
	bcftools call -c --ploidy 1  - | \
	vcfutils.pl vcf2fq | \
	seqtk seq -A | \
	bioawk -c fastx -v samp=$sample '{printf(">%s\n%s\n",samp,$seq)}'> ${sample}_consensus.fasta
}

export -f call_cons

parallel -j 40 call_cons 'NC_022830.1.fasta' ::: $(ls *_mitoreads.bam | tr '\n' ' ')

cat *_consensus.fasta > AllSamplesMitoConsensus.fasta