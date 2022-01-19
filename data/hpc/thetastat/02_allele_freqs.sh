module load parallel

# Since we don't know the ancestral state we use the ref 
# as the ancestral for Fst calculations as recommended
#
do_pop(){
	POP=$1

	REF=GCA_014634065.1_Adig_2.0_genomic.fna

    /fast/sci-irc/angsd/angsd -nThreads 32  -b ${POP}_bam.txt -anc $REF \
    -out ${POP}_af \
     -minMapQ 5 -minQ 20 -GL 2 -doSaf 1 \
     -rf chrs.txt \
     -sites sites.sorted.txt
}

export -f do_pop

parallel -j 3 do_pop ::: inshore north_offshore south_offshore

