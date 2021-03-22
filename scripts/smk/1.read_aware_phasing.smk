#!/usr/bin/env python

CHROMS=[i.strip() for i in open("scaffold_withSNPs.txt","r").readlines()]
extractPIRs="./bin/extractPIRs"
shapeit="./bin/shapeit"

rule all:
    input:
        expand("haplotypeData/{chr}.haps", chr=CHROMS),
        expand("haplotypeData/{chr}.sample", chr=CHROMS)
        #expand("PIRsList_by_chr/{chr}.PIRsList", chr=CHROMS)

rule split_by_scaffold:
    input:
        "Adigi.v2.filtered.indv74.vcf.gz"
    output:
        "vcf_by_chr/{chr}.vcf.gz"    
    shell:
        """
        vcftools --gzvcf {input} --max-missing 0.95 --chr {wildcards.chr} --recode --recode-INFO-all --stdout | bgzip > {output}
        """

rule extract_PIRs:
    input:
        vcf="vcf_by_chr/{chr}.vcf.gz",
        bamlist="bamlist_by_chr/{chr}"
    output:
        "PIRsList_by_chr/{chr}.PIRsList"
    shell:
        """
        {extractPIRs} --bam {input.bamlist} \
        --vcf {input.vcf} \
        --out {output} \
        --base-quality 20 \
        --read-quality 20
        """
"""
There will be some scaffolds can not be phased because of few samples without any genotype.
"""
rule phasing:
    input:
        vcf="vcf_by_chr/{chr}.vcf.gz",
        pir="PIRsList_by_chr/{chr}.PIRsList"
    output:
        haps="haplotypeData/{chr}.haps",
        sample="haplotypeData/{chr}.sample"
    shell:
        """
        {shapeit} -assemble \
        --input-vcf {input.vcf} \
        --input-pir {input.pir} \
        --thread 2 \
        --force \
        -O haplotypeData/{wildcards.chr}
        """
