align_uce(){
	in=$1
	mafft --adjustdirection --localpair --maxiterate 1000 $in | sed "s/>_R_/>/"> pomo_realignments/$(basename $in)
}

export -f align_uce

mkdir -p pomo_realignments

parallel -J 40 align_uce ::: pomo_alignments/*.fasta