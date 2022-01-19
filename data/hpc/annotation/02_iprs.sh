#conda activate interproscan

split_fasta(){
	input=$1
	cat $input | bioawk -c 'fastx' -v species_in=${input%.fasta} '{ if( (NR-1)%1000==0 ){file=sprintf("%s%d.fa",species_in,(NR-1));} printf(">%s\t%s\n%s\n",$name,$comment,$seq) >> file}'
}


f=protein.fasta
species=${f%.fasta}
echo ${species}
split_fasta $f
for seqs in ${species}*.fa;do
	/fast/sci-irc/interpro_scan/interproscan-5.53-87.0/interproscan.sh -i $seqs --disable-precalc -goterms --tempdir /fast/tmp/ > ${seqs%.fa}.log
done
