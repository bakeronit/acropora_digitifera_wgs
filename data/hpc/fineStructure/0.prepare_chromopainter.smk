#!/usr/bin/env python

CHROMS = [i.strip() for i in open('scaffolds_1M.txt','r').readlines()]

rule all:
    input:
        expand("recom/{chr}.recombfile",chr=CHROMS)

#step1: convert shapeit haps format to chromopainter
rule impute2chrmopainter:
    input:
        haps = "../phasing/haplotypeData/{chr}.haps"
    output:
        "phase_files/{chr}.phase" 
    shell:
        """
        perl ./fs_4.1.1/impute2chromopainter.pl {input.haps} {output}
        """
#step2: generate recombination map files
rule make_map:
    input:
        "phase_files/{chr}.phase"
    output:
        "recom/{chr}.recombfile"
    shell:
        """
        ./fs_4.1.1/makeuniformrecfile.pl {input} {output}
        """
