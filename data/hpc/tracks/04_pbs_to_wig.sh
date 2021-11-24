#!/bin/bash


pop="inshore"

cat fst.slidingwindow.tsv | awk -v pop=0 -f fst2wig.awk > ${pop}_pbs.wig
wigToBigWig ${pop}_pbs.wig adig2_chromsizes.txt ${pop}_pbs.bw



pop="north_offshore"

cat fst.slidingwindow.tsv | awk -v pop=1 -f fst2wig.awk > ${pop}_pbs.wig
wigToBigWig ${pop}_pbs.wig adig2_chromsizes.txt ${pop}_pbs.bw

pop="south_offshore"

cat fst.slidingwindow.tsv | awk -v pop=2 -f fst2wig.awk > ${pop}_pbs.wig
wigToBigWig ${pop}_pbs.wig adig2_chromsizes.txt ${pop}_pbs.bw
