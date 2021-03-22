#!/usr/bin/bash

window_pi(){
 pop=$1
 vcftools --gzvcf ../Adigi.v2.filtered.vcf.gz --keep ../${pop}.txt \
   --window-pi 10000 --window-pi-step 2000 --out ${pop}
}

export -f window_pi
parallel -j 3 window_pi ::: inshore northoffshore southoffshore

## remove windows with more than 30% bad sites

filter_windows(){
    pop=$1
    awk '{if($4<0.3) print $0}' ../0.window_missingness/${pop}.windows.missing.txt |\
    cut -f1,2,3 |grep -Fwf - ${pop}.windowed.pi > ${pop}.filtered_windowed.pi
}

export -f filter_windows
parallel -j 3 filter_windows ::: inshore northoffshore southoffshore
