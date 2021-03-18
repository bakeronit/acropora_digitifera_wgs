#!/usr/bin/bash

calculate_r2(){
    pop=$1
    shuf -n 20 ../${pop}.txt|awk '{print $0"\t"$0}' > ${pop}_subset.txt

    ../plink --vcf ../Adigi.v2.filtered.vcf.gz \
    --allow-no-sex --allow-extra-chr --double-id \
    --ld-window 999999 --ld-window-kb 100 --ld-window-r2 0 \
    --out ${pop} --thin 0.01 --r2 \
    --keep ${pop}_subset.txt
}

export -f calculate_r2

parallel -j 3 calculate_r2 ::: inshore northoffshore southoffshore


