#!/bin/bash

for f in *.pestPG;do
	pop=${f%.pestPG}
	cat $f | awk -f theta2wig.awk > ${pop}_TajimaD.wig

	wigToBigWig ${pop}_TajimaD.wig adig2_chromsizes.txt ${pop}_TajimaD.bw
done

