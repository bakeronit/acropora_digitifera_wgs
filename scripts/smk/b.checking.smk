#!/usr/bin/bash

PREFIX="Adigi.indel5bp_snp_hardfilter_passed_biallelic_mdust"
outprefix="5.indel5bp_snp_hardfilter_passed_biallelic_mdust"

rule all:
    input:
        outprefix + "_depth.idepth",
        outprefix + "_site_depth.ldepth",
        outprefix + "_site_depth_mean.ldepth.mean",
        outprefix + "_missing_indv.imiss",
        outprefix + "_missing_site.lmiss",
        outprefix + "_relatedness.relatedness2",

"""
Generate a file containing the mean depth per individual (suffixed '.idepth')
"""

rule depth:
    input:
        PREFIX + ".vcf.gz"
    output:
        outfile=outprefix + "_depth.idepth",
    shell:
        """
        vcftools --gzvcf {input} --depth --out {outprefix}_depth  
        """

rule site_depth:
    input:
        PREFIX + ".vcf.gz"
    output:
        outfile = outprefix + "_site_depth.ldepth",
    shell:
        """
        vcftools --gzvcf {input} --site-depth --out {outprefix}_site_depth
        """

rule site_depth_mean:
    input:
        PREFIX + ".vcf.gz"
    output:
        outfile = outprefix + "_site_depth_mean.ldepth.mean",
    shell:
        """
        vcftools --gzvcf {input} --site-mean-depth --out {outprefix}_site_depth_mean
        """

rule missing_indv:
    input:
        PREFIX + ".vcf.gz"
    output:
        outfile = outprefix + "_missing_indv.imiss",
    shell:
        """
        vcftools --gzvcf {input} --missing-indv --out {outprefix}_missing_indv
        """

rule missing_site:
    input:
        PREFIX + ".vcf.gz"
    output:
        outfile = outprefix + "_missing_site.lmiss",
    shell:
        """
        vcftools --gzvcf {input} --missing-site --out {outprefix}_missing_site
        """

rule relatedness:
    input:
        PREFIX + ".vcf.gz"
    output:
        outfile = outprefix + "_relatedness.relatedness2",
    shell:
        """
        vcftools --gzvcf {input} --relatedness2 --out {outprefix}_relatedness
        """
