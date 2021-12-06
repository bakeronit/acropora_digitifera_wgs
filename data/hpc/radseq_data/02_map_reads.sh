#!/bin/bash

#conda activate gatk


GENOME=GCA_014634065.1_Adig_2.0_genomic.fna

export PATH=${PATH}:/fast/sci-irc/abalone_relatedness/hpc/gatk/gatk-4.2.1.0/

#set -e

# See here 
# http://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups#latest
# For a discussion of assigning read groups etc
#
fastq2ubam(){
	input1=$1
	output=$2
	sample=$3
	flowcell=$4
	lane=$5

	gatk FastqToSam \
		-FASTQ $input1 \
    	-OUTPUT $output.bam \
    	-READ_GROUP_NAME $sample.$flowcell.$lane \
    	-SAMPLE_NAME $sample \
    	-LIBRARY_NAME $sample \
    	-PLATFORM_UNIT $flowcell.$lane.$sample \
    	-PLATFORM ILLUMINA
}

markadapters(){
	input=$1
	output=$2

    gatk MarkIlluminaAdapters \
    -I $input.bam \
    -M ${input}_txt \
    -O $output.bam 
}

map_reads(){
	input=$1
	fasta=$2
	ubam=$3
	output=$4

	gatk SamToFastq \
    -I $input.bam \
    -FASTQ /dev/stdout \
    -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
    |  \
    bwa mem -M -t 16 -p $fasta /dev/stdin \
    | \
    gatk MergeBamAlignment \
    -ALIGNED_BAM /dev/stdin \
    -UNMAPPED_BAM $ubam \
    -OUTPUT $output.bam \
    -R $fasta -CREATE_INDEX true -ADD_MATE_CIGAR true \
    -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true \
    -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 \
    -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS 
}

# First index the genome
#bwa index $GENOME

# Make the dict
#gatk CreateSequenceDictionary -REFERENCE $GENOME -OUTPUT ${GENOME%.fna}.dict


for f in *.FASTQ.gz; do

	flowcell=$(gunzip -c $f | head -n 1 | awk -F ':' '{print $3}')
	lane=$(gunzip -c $f | head -n 1 | awk -F ':' '{print $4}')
	sample=${f%.FASTQ.gz}

	input1=$f

	echo $input1 $sample $sample $flowcell $lane

	if [ ! -f $sample.bam ]; then
		fastq2ubam $input1 ${sample} $sample $flowcell $lane
	fi

	if [ ! -f ${sample}_markadapters.bam ]; then
		markadapters $sample ${sample}_markadapters
	fi

	if [ ! -f ${sample}_mapped.bam ]; then
		map_reads ${sample}_markadapters $GENOME $sample.bam ${sample}_mapped
		samtools index ${sample}_mapped.bam
	fi
done




