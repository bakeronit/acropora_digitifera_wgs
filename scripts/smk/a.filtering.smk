#!/usr/bin/python

REFERENCE="genome/reference.fa"

rule all:
    input:
        "0.gatk_raw.stats",
        "1.indel5bp.stats",
        "2.indel5bp_snp.stats",
        "3.indel5bp_snp_hardfilter.stats",
        "4.indel5bp_snp_hardfilter_passed_biallelic.stats",
        "5.indel5bp_snp_hardfilter_passed_biallelic_mdust.stats",
        "Adigi.indel5bp_snp_hardfilter_passed_biallelic_mdust.vcf.gz"


rule indel5bp:
    input:
        "Adigi.gatk_raw.vcf.gz"
    output:
        vcf=temp("Adigi.indel5bp.vcf.gz"),
        control="0.gatk_raw.stats",
        stats="1.indel5bp.stats"
    threads: 20
    shell:
        """
        bcftools stats --threads {threads} {input} > {output.control}
        bcftools filter -g 5 --threads {threads} -O z -o {output.vcf} {input}
        tabix -p vcf {output.vcf}
        
        bcftools stats --threads {threads} {output.vcf} > {output.stats}
        """

rule gatkSelectSNPs:
    input:
        "Adigi.indel5bp.vcf.gz"
    output:
        vcf=temp("Adigi.indel5bp_snp.vcf.gz"),
        stats="2.indel5bp_snp.stats"
    threads: 20
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type SNP \
        -O {output.vcf}

        bcftools stats --threads {threads} {output.vcf} > {output.stats}
        """

rule gatkHardfilter:
    input:
        "Adigi.indel5bp_snp.vcf.gz"
    output:
        vcf=temp("Adigi.indel5bp_snp_hardfilter.vcf.gz"),
        stats="3.indel5bp_snp_hardfilter.stats"
    threads: 20
    shell:
        """
        gatk VariantFiltration \
        -V {input} \
        -filter "QD < 10.0" --filter-name "QD10" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O {output.vcf}

        bcftools stats --threads {threads} {output.vcf} > {output.stats}
        """

rule gatkSelectPassedSNPs:
    input:
        vcf="Adigi.indel5bp_snp_hardfilter.vcf.gz"
    output:
        vcf=temp("Adigi.indel5bp_snp_hardfilter_passed_biallelic.vcf.gz"),
        stats="4.indel5bp_snp_hardfilter_passed_biallelic.stats"
    threads: 20
    shell:
        """
        gatk SelectVariants \
        -V {input.vcf} \
        -O {output.vcf} \
        --restrict-alleles-to BIALLELIC \
        --exclude-filtered true

        bcftools stats --threads {threads} {output.vcf} > {output.stats}
        """

rule removeSimpleRepeatRegion:
    input:
        vcf="Adigi.indel5bp_snp_hardfilter_passed_biallelic.vcf.gz",
        reference=REFERENCE
    output:
        bed="genome.mdust.bed",
        vcf="Adigi.indel5bp_snp_hardfilter_passed_biallelic_mdust.vcf.gz",
        stats="5.indel5bp_snp_hardfilter_passed_biallelic_mdust.stats"
    threads: 20
    shell:
        """
        ./bin/mdust {input.reference} -c |cut -f1,3,4 > {output.bed}

        vcftools --gzvcf {input.vcf} \
        --exclude-bed {output.bed} \
        --recode --recode-INFO-all --stdout | \
        bgzip > {output.vcf}

        bcftools stats --threads {threads} {output.vcf} > {output.stats}
        """
