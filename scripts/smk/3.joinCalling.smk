#!/usr/bin/env python

"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""

REFERENCE="data/genome/reference.fa"
scaffolds = []
SAMPLES=[i.strip() for i in open('samples.txt','r').readlines()]

with open(REFERENCE, 'rt') as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            line = line.split(" ")[0]
            scaffolds.append(line[1:])

wildcard_constraints:
    sample = "|".join(SAMPLES),
    scaffold = "|".join(scaffolds)

rule all:
    input:
        "analysis/vcf/all_raw.vcf.gz",
        "analysis/vcf/all_raw.vcf.gz.tbi"

rule GenomicsDBImport:
    input:
        "sample_map.txt"
    output:
        db = directory("analysis/genomicsDB/{scaffold}.db"),
        tar = "analysis/genomicsDB/{scaffold}.db.tar"
    threads: 6
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output.db} \
        --L {wildcards.scaffold} \
        --sample-name-map {input} \
        --reader-threads 5 \
        --batch-size 50

        tar -cf {output.tar} {output.db}
        """

## --use-new-qual-calculator \  # deprecated after gatk 4.1.4
rule GenotypeGVCFs:
    input:
        "analysis/genomicsDB/{scaffold}.db"
    output:
        "analysis/vcf/{scaffold}.vcf.gz"
    shell:
        """
        gatk --java-options "-Xmx8g -Xms8g" \
        GenotypeGVCFs \
        -R {REFERENCE} \
        -O {output} \
        --only-output-calls-starting-in-intervals \
        -V gendb://{input} \
        -L {wildcards.scaffold}
        """

rule GatherVcfs:
    input:
        expand("analysis/vcf/{scaffold}.vcf.gz",scaffold = scaffolds)
    output:
        file="analysis/vcf/all_raw.vcf.gz",
        index="analysis/vcf/all_raw.vcf.gz.tbi"
    params:
        " -I ".join("analysis/vcf/" + s + ".vcf.gz" for s in scaffolds)
    shell:
        """
        gatk --java-options "-Xmx6g -Xms6g" \
        GatherVcfsCloud \
        --ignore-safety-checks \
        --gather-type BLOCK \
        -I {params} \
        -O {output.file}

        gatk --java-options "-Xmx6g -Xms6g" \
        IndexFeatureFile \
        -I {output.file} \
        -O {output.index}
        """

