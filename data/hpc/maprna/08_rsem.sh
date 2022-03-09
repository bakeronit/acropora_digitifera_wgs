
# rsem-prepare-reference --gff3 ../../genome/adig-v2-ncbi.gff \
#                             --star \
#                             -p 8 \
#                             GCA_014634065.1_Adig_2.0_genomic.fna starref
                            

do_rsem(){
    f1=$1
    sample=${f1%_1.fastq}
    echo $sample
    rsem-calculate-expression --paired-end --alignments -p 8 \
        ${sample}Aligned.toTranscriptome.out.bam ./starref ${sample}        
}

export -f do_rsem

parallel -j 12 do_rsem ::: *_1.fastq
