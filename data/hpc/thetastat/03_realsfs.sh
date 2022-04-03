module load parallel

# Since we don't know the ancestral state we use the ref as the ancestral for Fst calculations
#
do_realsfs(){
	POP=$1

	/fast/sci-irc/angsd/misc/realSFS ${POP}_af.saf.idx -P 24 -fold 1 > ${POP}.sfs
	/fast/sci-irc/angsd/misc/realSFS saf2theta ${POP}_af.saf.idx -outname ${POP} -sfs ${POP}.sfs -fold 1

}

export -f do_realsfs

do_realsfs jp

#parallel -j 3 do_realsfs ::: inshore north_offshore south_offshore

