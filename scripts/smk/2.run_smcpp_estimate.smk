#!/usr/bin/env python

POPS=["inshore","northoffshore","southoffshore"]
CHROMS=[i.strip() for i in open("scaffolds_1M.txt",'r').readlines()]

rule all:
    input:
        expand("estimate/{pop}.final.json", pop=POPS)

rule smc:
    input:
        expand("smc/{{pop}}/{chr}.smc.gz",chr=CHROMS)
    output:
        "estimate/{pop}.final.json"
    threads: 30
    shell:
        """
        smc++ estimate --cores {threads} -o estimate --base {wildcards.pop} \
        --timepoints 200 200000 --em-iterations 50 --thinning 2000 --knots 40 1.2e-8 {input} 
        """

