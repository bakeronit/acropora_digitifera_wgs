#!/usr/bin/env python

SAMPLES=[i.strip() for i in open('samples.txt','r').readlines()]

rule all:
    input:
        expand("data/ubam/{sample}.unaligned_reads.bam", sample=SAMPLES)

rule fastq2ubam:
    input:
        r1 = "data/fastq/{sample}_R1_001.fastq.gz",
        r2 = "data/fastq/{sample}_R2_001.fastq.gz"
    output:
        ubam = "data/ubam/{sample}.unaligned_reads.bam"
    threads: 1
    shell:
        """
        picard FastqToSam \
          FASTQ={input.r1} \
          FASTQ2={input.r2} \
          USE_SEQUENTIAL_FASTQS=true \
          OUTPUT={output.ubam} \
          PLATFORM=illumina \
          SAMPLE_NAME={wildcards.sample} \
          PLATFORM_MODEL=NovaSeq6000_S4
        """
