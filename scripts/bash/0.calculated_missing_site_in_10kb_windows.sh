#!/usr/bin/bash

path=../../01.mapping_genotyping_variantcalling/analysis/mapping
declare -A BAM
BAM[inshore]=BR_4_091_S113_L004_aligned_duplicates_marked_sorted.bam
BAM[northoffshore]=AR_125_388_S128_L004_aligned_duplicates_marked_sorted.bam
BAM[southoffshore]=RS2_C11_784_S156_L004_aligned_duplicates_marked_sorted.bam

bioawk -c fastx '{print $name"\t1\t"length($seq)}' \
 ../../00.rawdata/genome/reference.fa > reference.bed

#create windows
bedtools makewindows -b reference.bed -w 9999 -s 2000 > windows_10k.bed

# calculate window called site missingness
cal_missing(){
    pop=$1  
    bedtools genomecov -ibam ${path}/${BAM[$pop]} -d | \
    awk '{if($3<3) print$0}' | \
    awk '{print $1"\t"$2"\t"$2}' |\
    bedtools merge |\
    bedtools intersect -a windows_10k.bed -b - -wao | \
    cut -f1,2,3,7 |\
    awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' |\
    sort -k1,1 -k2,2n -k3,3n > ${pop}.windows.missing.txt
}

export -f cal_missing
cal_missing inshore 
cal_missing northoffshore 
cal_missing southoffshore 

