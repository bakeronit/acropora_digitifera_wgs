module load mafft/7.394
module load samtools/1.11
module load parallel/20180122

# This creates common_uces.bed
# R CMD BATCH find_common_uces.R

rm alignments/*


gather_uce(){
	cuce=$1
	echo "Common UCE $cuce"



	adi_uces=$(cat common_uces.tsv | awk -v uce=$cuce '($4==uce) && ($6=="adi"){print $5}' | sort -u | tr '\n' ' ')
	aten_uces=$(cat common_uces.tsv | awk -v uce=$cuce '($4==uce) && ($6=="aten"){print $5}' | sort -u | tr '\n' ' ')
	amil_uces=$(cat common_uces.tsv | awk -v uce=$cuce '($4==uce) && ($6=="amil"){print $5}' | sort -u | tr '\n' ' ')		
	echo "${adi_uces} ${aten_uces} ${amil_uces}"

	samtools faidx adi.uces.fasta ${adi_uces} | sed "s/>.*/>adi/" >> alignments/${cuce}.fasta

	for ucefasta in adi_*.uces.fasta;do
		id=${ucefasta%.uces.fasta}
		samtools faidx $ucefasta ${adi_uces} | sed "s/>.*/>${id}/" >> alignments/${cuce}.fasta
	done

	for ucefasta in [FM]I*.uces.fasta;do
		id=${ucefasta%.uces.fasta}
		samtools faidx $ucefasta ${aten_uces} | sed "s/>.*/>${id}/" >> alignments/${cuce}.fasta
	done


	# Amil
	samtools faidx amil.uces.fasta ${amil_uces} | sed "s/>.*/>amil/" >> alignments/${cuce}.fasta

}

align_uce(){
	cuce=$1
	mafft --adjustdirection --localpair --maxiterate 1000 alignments/${cuce}.fasta | sed "s/>_R_/>/"> alignments/${cuce}_align.fasta
}

export -f gather_uce
export -f align_uce

parallel -J 40 gather_uce ::: $(cat common_uces.tsv | awk '{print $4}' | sort -u | tr '\n' ' ')

parallel -J 40 align_uce ::: $(cat common_uces.tsv | awk '{print $4}' | sort -u | tr '\n' ' ')


