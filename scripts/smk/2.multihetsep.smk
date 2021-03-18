#!/usr/bin/env python

SAMPLES = [ i.strip() for i in open("samples_12.txt","r").readlines() ]
#CHROMS=[ i.strip() for i in open("scaffolds_1M.txt","r").readlines() ][2:5]
CHROMS=['BLFC01000927.1']
rule all:
    input:
        expand("multihetsep/indv12.{chr}.multihetsep.txt",chr=CHROMS)

rule generate:
    input:
        mask = expand("samples/{sample}/{{chr}}_mask.bed.gz",sample=SAMPLES),
        vcf = expand("samples/{sample}/{{chr}}.vcf.gz",sample=SAMPLES)
    output:
        "multihetsep/indv12.{chr}.multihetsep.txt"
    shell:
        """
        ./bin/msmc-tools/generate_multihetsep.py --chr {wildcards.chr} \
        --mask {input.mask[0]} \
        --mask {input.mask[1]} \
        --mask {input.mask[2]} \
        --mask {input.mask[3]} \
        --mask {input.mask[4]} \
        --mask {input.mask[5]} \
        --mask {input.mask[6]} \
        --mask {input.mask[7]} \
        --mask {input.mask[8]} \
        --mask {input.mask[9]} \
        --mask {input.mask[10]} \
        --mask {input.mask[11]} \
        {input.vcf[0]}  \
        {input.vcf[1]}  \
        {input.vcf[2]}  \
        {input.vcf[3]}  \
        {input.vcf[4]}  \
        {input.vcf[5]}  \
        {input.vcf[6]}  \
        {input.vcf[7]}  \
        {input.vcf[8]}  \
        {input.vcf[9]}  \
        {input.vcf[10]}  \
        {input.vcf[11]}  > {output}
        """
