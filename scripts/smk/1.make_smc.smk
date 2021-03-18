#!/usr/bin/env python

CHROMS=[i.strip() for i in open("scaffolds_1M.txt",'r').readlines()]
POPS = ["inshore","northoffshore","southoffshore"]

rule all:
    input:
        expand("smc/{pop}/{chr}.smc.gz", pop=POPS, chr=CHROMS)


rule bcf_query:
    input:
        "Adigi.v2.filtered.vcf.gz"
    output:
        expand("{pop}.txt",pop=POPS)
    run:
        inshore = output[0]
        noffshore = output[1]
        soffshore = output[2]
        shell("bcftools query -l {input} |grep -E 'AI|BR' > {inshore}")
        shell("bcftools query -l {input} |grep 'AR' > {noffshore}")
        shell("bcftools query -l {input} |grep 'RS' > {soffshore}")

rule split_vcf:
    input:
        "Adigi.v2.filtered.vcf.gz"
    output:
        vcf = "vcf_files/{chr}.vcf.gz",
        tbi = "vcf_files/{chr}.vcf.gz.tbi"
    shell:
        """
        bcftools view -r {wildcards.chr} -Oz -o {output.vcf} {input}
        tabix {output.vcf}
        """

rule convert:
    input:
        vcf = "vcf_files/{chr}.vcf.gz",
        sample = "{pop}.txt"
    output:
        "smc/{pop}/{chr}.smc.gz"
    params: 
        prefix = lambda wildcards: wildcards.pop[:2].upper()
    shell:
        """
        smc++ vcf2smc {input.vcf} {output} {wildcards.chr} \
        {params.prefix}:$(paste -s -d ',' {input.sample} )
        """
