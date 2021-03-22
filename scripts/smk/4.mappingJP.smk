#!/usr/bin/env python

"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""
SAMPLES=[i.strip() for i in open('JP_samples.txt','r').readlines()]
REFERENCE="data/genome/reference.fa"

rule all:
    input:
        expand("analysis/chuya_mapping/{sample}_aligned_duplicates_marked_sorted.bam",sample=SAMPLES)

rule fastq2ubam:
    input:
        r1="data/chuya_fastq/{sample}_1.fastq.gz",
        r2="data/chuya_fastq/{sample}_2.fastq.gz"
    output:
        "data/chuya_ubam/{sample}.unaligned_reads.bam"
    shell:
        """
        picard FastqToSam \
          FASTQ={input.r1} \
          FASTQ2={input.r2} \
          OUTPUT={output} \
          PLATFORM=illumina \
          SAMPLE_NAME={wildcards.sample}
        """

rule SamToFastqAndBwaMem:
    input:
        ubam = "data/chuya_ubam/{sample}.unaligned_reads.bam",
        reference=REFERENCE
    output:
        "analysis/chuya_mapping/{sample}_unmerged.bam"
    log:
        "logs/{sample}.bwa.stderr.log"
    threads: 15
    shell:
        """
        java -Dsamjdk.compression_level=5 -Xms3000m -jar $PICARD_HOME/picard.jar SamToFastq \
        INPUT={input.ubam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
        bwa mem -K 100000000 -p -v 3 -t {threads} -Y {input.reference} /dev/stdin - 2> >(tee {log} >&2) | \
        samtools view -1 - > {output}
        """

bwa_version = subprocess.check_output("bwa 2>&1 | grep -e 'Version'", shell=True).decode("utf-8").rstrip()

rule MergeBamAlignment:
    input:
        unmerged_bam = "analysis/chuya_mapping/{sample}_unmerged.bam",
        unmapped_bam = "data/chuya_ubam/{sample}.unaligned_reads.bam",
        reference = REFERENCE
    output:
        temp("/scratch/jc502059/chuya_mapping/{sample}_aligned_unsorted.bam")
    threads: 1 
    params:
        v_bwa = bwa_version

    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ALIGNED_BAM {input.unmerged_bam}  \
        --UNMAPPED_BAM {input.unmapped_bam} \
        --OUTPUT {output} \
        --REFERENCE_SEQUENCE {input.reference} \
        --PAIRED_RUN true \
        --SORT_ORDER "unsorted" \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --PROGRAM_RECORD_ID "bwamem" \
        --PROGRAM_GROUP_VERSION "{params.v_bwa}" \
        --PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t 2 -Y {input.reference}" \
        --PROGRAM_GROUP_NAME "bwamem" \
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --UNMAP_CONTAMINANT_READS true
        """

rule MarkDuplicates:
    input:
        "/scratch/jc502059/chuya_mapping/{sample}_aligned_unsorted.bam"
    output:
        bam = temp("/scratch/jc502059/chuya_mapping/{sample}_aligned_unsorted_duplicates_marked.bam"),
        metrics_filename = "analysis/chuya_mapping/{sample}.duplicate_metrics"
    threads: 1
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m -XX:+UseParallelGC -XX:ParallelGCThreads=2" \
        MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics_filename} \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER "queryname" \
        --CREATE_MD5_FILE true
        """

rule SortAndFixTags:
    input:
        bam = "/scratch/jc502059/chuya_mapping/{sample}_aligned_unsorted_duplicates_marked.bam",
        reference = REFERENCE
    output:
        protected("analysis/chuya_mapping/{sample}_aligned_duplicates_marked_sorted.bam"),
    threads: 1
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" \
        SortSam \
        --INPUT {input.bam} \
        --OUTPUT /dev/stdout \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX false \
        --CREATE_MD5_FILE false \
        | \
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" \
        SetNmMdAndUqTags \
        --INPUT /dev/stdin \
        --OUTPUT {output} \
        --CREATE_INDEX true \
        --CREATE_MD5_FILE true \
        --REFERENCE_SEQUENCE {input.reference} 
        """
