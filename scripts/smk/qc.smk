#!/usr/bin/env python

SAMPLES=[i.strip() for i in open('samples.txt','r').readlines()]
rule all:
    input:
        expand("analysis/qc/bamtools/{sample}_bamtools.stats",sample=SAMPLES),
        expand("analysis/qc/coverage/{sample}_coverage.txt",sample=SAMPLES)
rule bamstats:
    input:
        "analysis/mapping/{sample}_aligned_duplicates_marked_sorted.bam"
    output:
        "analysis/qc/bamtools/{sample}_bamtools.stats"
    shell:
        """
        module load bamtools
        bamtools stats -in {input} | grep -v "*" > {output}
        """

rule genomeCov:
    input:
        "analysis/mapping/{sample}_aligned_duplicates_marked_sorted.bam"
    output:
        "analysis/qc/coverage/{sample}_coverage.txt"
    shell:
        """
        bamcov {input} -o {output}
        """
