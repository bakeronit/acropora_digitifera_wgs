#!/usr/bin/bash

cal_td(){
    pop=$1
    bcftools view -S ../${pop}.txt ../Adigi.v2.filtered.vcf.gz |\
    vk tajima 10000 2000 - |\
    sed '1d' |\
    awk '{print $1"\t"$2+1"\t"$3"\t"$6}' > ${pop}_1based_td.txt
}

#export -f cal_td
#parallel -j 3 cal_td ::: inshore northoffshore southoffshore

filter_td(){
    pop=$1
    awk '{if($4<0.3)print $0}' ../0.window_missingness/${pop}.windows.missing.txt |\
    cut -f1,2,3 |\
    grep -Fwf - ${pop}_1based_td.txt > ${pop}.filtered_windowed.td
}

export -f filter_td
parallel -j 3 filter_td ::: inshore northoffshore southoffshore
