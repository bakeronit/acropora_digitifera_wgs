#!/usr/bin/bash

vcf=../Adigi.v2.filtered.vcf.gz
cal_fst(){
    pop1=$1
    pop2=$2
    vcftools --gzvcf $vcf \
    --weir-fst-pop ../${pop1}.txt \
    --weir-fst-pop ../${pop2}.txt \
    --out ${pop1}_${pop2}_per_snp_fst
}

export -f cal_fst

cal_fst inshore northoffshore
cal_fst inshore southoffshore
cal_fst northoffshore southoffshore
