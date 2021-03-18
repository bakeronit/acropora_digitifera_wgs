#!/usr/bin/env python

CHROMS=[i.strip() for i in open("scaffolds_1M.txt",'r').readlines()]
rule all:
    input:
        expand("split/{pop1}_{pop2}.final.json",zip,pop1=["inshore","inshore","northoffshore"],pop2=["northoffshore","southoffshore","southoffshore"])

rule vcf2smc_jfs:
    input:
        vcf="vcf_files/{chr}.vcf.gz",
        pop1="{pop1}.txt",
        pop2="{pop2}.txt"
    output:
        "smc/jfs_{pop1}_{pop2}/{chr}.smc.gz"
    params:
        p1=lambda wildcards: wildcards.pop1[:2].upper(),
        p2=lambda wildcards: wildcards.pop2[:2].upper()
    shell:
        """
        smc++ vcf2smc {input.vcf} {output} {wildcards.chr} \
        {params.p1}:$(paste -s -d ',' {input.pop1}) \
        {params.p2}:$(paste -s -d ',' {input.pop2})
        """

rule split:
    input:
        estimate1="estimate/{pop1}.final.json",
        estimate2="estimate/{pop2}.final.json",
        jfs = expand("smc/jfs_{{pop1}}_{{pop2}}/{chr}.smc.gz",chr=CHROMS)
    output:
        "split/{pop1}_{pop2}.final.json"
    threads: 30
    shell:
        """
        smc++ split -o split --base {wildcards.pop1}_{wildcards.pop2} \
        --timepoints 200 200000 --em-iterations 50 --thinning 2000 \
        {input.estimate1} {input.estimate2} {input.jfs}
        """
