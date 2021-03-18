#!/usr/bin/env python

CHROMS=[i.strip() for i in open("scaffold_phased.txt","r").readlines()]
shapeit="./bin/shapeit"

rule all:
    input:
        expand("phased_vcf/{chr}.vcf.gz",chr=CHROMS)

rule convert:
    input:
        haps="haplotypeData/{chr}.haps",
        sample="haplotypeData/{chr}.sample"
    output:
        "phased_vcf/{chr}.vcf.gz"
    shell:
        """
        {shapeit} -convert \
        --input-haps haplotypeData/{wildcards.chr} \
        --output-vcf phased_vcf/{wildcards.chr}.vcf

        bgzip phased_vcf/{wildcards.chr}.vcf
        tabix {output}
        """

