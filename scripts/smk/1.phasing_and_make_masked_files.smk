#!/usr/bin/env python3

SAMPLES = [ i.strip() for i in open("samples_12.txt","r").readlines() ]
SCAFFOLDS = [ i.strip() for i in open("scaffolds_1M.txt","r").readlines() ]
REFERENCE = "reference.fa"
MSMC_HOME = "/home/5/jc502059/bioprojects/adigitifera_wgs/05.MSMC_phased_new/bin"
EXTRACT_PIRs = "/home/5/jc502059/bioprojects/adigitifera_wgs/bin/extractPIRs.v1.r68.x86_64/extractPIRs"
SHAPEIT = "/home/5/jc502059/bioprojects/adigitifera_wgs/bin/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"

rule all:
    input:
        expand("samples/{sample}/{scaffold}_mask.bed.gz", sample = SAMPLES,scaffold = SCAFFOLDS ),
        expand("samples/{sample}/{scaffold}.vcf.gz", sample = SAMPLES, scaffold = SCAFFOLDS),
        #expand("finalVCF/{scaffold}.vcf.gz", scaffold = SCAFFOLDS),
        #expand("mappability/Adigi_{scaffold}.mask.bed.gz",scaffold = SCAFFOLDS)

rule call:
    input:
        reference = REFERENCE,
        bamfile = "bamfiles.txt"
    output:
        "vcf_files/{scaffold}.vcf.gz"
    shell:
        "bcftools mpileup -q 30 -Q 30 -C 50 -Oz -r {wildcards.scaffold} -f {input.reference} $(cat {input.bamfile} | tr '\n' ' ') |\
        bcftools call -c -V indels | bcftools view -M 2 -Oz > {output}"

rule extract_PIRs:
    input:
        bamlist = "bamlist_files/{scaffold}",
        vcf = "vcf_files/{scaffold}.vcf.gz"
    output:
        "extractPIRs/{scaffold}.PIRsList"
    shell:
        """
        {EXTRACT_PIRs} --bam {input.bamlist} \
        --vcf {input.vcf} \
        --out {output} \
        --base-quality 20 --read-quality 20 
        """

rule assemble:
    input:
        vcf = "vcf_files/{scaffold}.vcf.gz",
        pir = "extractPIRs/{scaffold}.PIRsList"
    output:
        "haplotypeData/{scaffold}.haps",
    params:
        prefix=lambda wildcard, output: output[0][:-5]
    threads: 20
    shell:
        """
        {SHAPEIT} -assemble \
        --thread {threads} \
        --input-vcf {input.vcf} \
        --input-pir {input.pir} \
        -O {params.prefix} 
        """

rule convert:
    input:
        "haplotypeData/{scaffold}.haps"
    output:
        phased_vcf = "phased_VCF/{scaffold}.phased.vcf"
    params:
        prefix=lambda wildcard, input: input[0][:-5]
    shell:
        """
        {SHAPEIT} -convert \
        --input-haps {params.prefix} \
        --output-vcf {output.phased_vcf}
        """

rule zip:
    input:
        phased = "phased_VCF/{scaffold}.phased.vcf",
        vcf = "vcf_files/{scaffold}.vcf.gz"
    output:
        "phased_VCF/{scaffold}.phased.vcf.gz"
    shell:
        """
        bcftools view -Oz {input.phased} > {output}
        bcftools index -f {output}
        bcftools index -f {input.vcf} 
        """

rule merge:
    input:
        phased = "phased_VCF/{scaffold}.phased.vcf.gz",
        vcf = "vcf_files/{scaffold}.vcf.gz"
    output:
        "finalVCF/{scaffold}.vcf.gz"
    shell:
        """
        bcftools merge --force-samples {input.vcf} {input.phased} |\
        awk 'BEGIN {{OFS="\t"}}
            $0 ~/^##/ {{print $0}}
            $0 ~/^#CHROM/ {{for(i=1;i<84;i++) printf "%s"OFS, $i; print $84}}
            $0 !~/^#/ {{ for(x=10;x<=84;x++) $x=$(x+75); for(j=1;j<84;j++) printf "%s"OFS, $j; print $84 }}' | \
            bcftools view -Oz > {output}
            tabix -p vcf {output}
        """

rule write_newheader:
    input:
        vcf = "finalVCF/{scaffold}.vcf.gz"
    output:
        head = temp("{scaffold}.head")
    run:
        import subprocess
        out = open(output[0],'w')
        for line in shell("bcftools view -h {input.vcf}", iterable=True):
            line = line.strip()
            if line.startswith("##"):
                print(line, file=out)
            else:
                samples = []
                new_line = line.split()[:9]
                for id in line.split()[9:]:
                    newid = id.split("_")[:3]
                    samples.append("_".join(newid))
                new_line.extend(samples)
                print("\t".join(new_line), file=out)

rule reheader:
    input:
        head = "{scaffold}.head",
        vcf = "finalVCF/{scaffold}.vcf.gz"
    output:
        "finalVCF/{scaffold}.reheader.vcf.gz"
    shell:
        """
        bcftools reheader -h {input.head} {input.vcf} > {output}
        """


rule bamcaller:
    input:
        reference = REFERENCE,
        bam = "bam/{sample}.bam",
        vcf = "finalVCF/{scaffold}.reheader.vcf.gz"
    output:
        bed = "samples/{sample}/{scaffold}_mask.bed.gz",
        vcf = "samples/{sample}/{scaffold}.vcf.gz"

    shell:
        """
        mean_cov=$(samtools depth -r {wildcards.scaffold} {input.bam}| awk '{{sum += $3}} END {{print sum/NR}}')
        
        bcftools view -s {wildcards.sample} {input.vcf} | \
        {MSMC_HOME}/msmc-tools/bamCaller.py $mean_cov {output.bed} | gzip -c > {output.vcf}
        """

rule mask_fa:
    input:
        reference = REFERENCE
    output:
        raw_fasta = "mappability/Adigi.rawMask_35.fa",
        fasta = "mappability/Adigi.mask_35_50.fa"
    shell:
        """
        ./bin/seqbility-20091110/splitfa {input.reference} 35 | split -l 20000000
        ls x* | while read id; do bwa aln -t 20 -R 1000000 -O 3 -E 3 {input.reference} $id > ${{id}}.sai;done
        ls x??| while read id; do bwa samse -f ${{id}}.sam {input.reference} ${{id}}.sai ${{id}} && gzip ${{id}}.sam;done
        wait
        gzip -dc x??.sam.gz | ./bin/seqbility-20091110/gen_raw_mask.pl > {output.raw_fasta} 
        ./bin/seqbility-20091110/gen_mask -l 35 -r 0.5 {output.raw_fasta} > {output.fasta}
        """

rule mappability:
    input:
        "mappability/Adigi.mask_35_50.fa"
    output:
        expand("mappability/Adigi_{scaffold}.mask.bed.gz",scaffold = SCAFFOLDS)
    shell:
        """
        python2 ./bin/msmc-tools/makeMappabilityMask.py
        """

