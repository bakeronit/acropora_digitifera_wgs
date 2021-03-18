#!/usr/bin/env python


CHROMS = [i.strip() for i in open("scaffold_ids.txt","r").readlines()]
selscan = "./bin/selscan"

rule all:
    input:
        #expand("vcf_files/{pop}/{chr}.vcf.gz", pop = ["inshore","northoffshore","southoffshore"],chr=CHROMS),
        #expand("map_files/{chr}.map.txt",chr=CHROMS)
        expand("ihs_out/{pop}/{chr}.ihs.out",pop=["inshore","northoffshore","southoffshore"],chr=CHROMS)


rule extract:
    input:
        "Adigi.v2.indv74_phased.vcf.gz"
    output:
        "vcf_files/{pop}/{chr}.vcf.gz"
    shell:
        """
        bcftools view -S {wildcards.pop}.txt -r {wildcards.chr} {input} |gzip > {output}
        """

rule make_map:
    input:
        "../03.phasing/phased_vcf/{chr}.vcf.gz"
    output:
        "map_files/{chr}.map.txt"
    run:
        import gzip
        cM = 3.6
        out = open(output[0],'w')
        with gzip.open(input[0],'rt') as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                sid, pos = line.strip().split()[:2]
                gpos = cM/1000000*int(pos)
                print(f'{sid}\t{sid}:{pos}\t{gpos}\t{pos}', file = out)

rule selscan_ihs:
    input:
        vcf = "vcf_files/{pop}/{chr}.vcf.gz",
        map = "map_files/{chr}.map.txt"
    output:
        "ihs_out/{pop}/{chr}.ihs.out"
    threads: 10
    shell:
        """
        {selscan} --ihs --vcf {input.vcf} --map {input.map} \
        --threads {threads} --out ihs_out/{wildcards.pop}/{wildcards.chr}
        """

