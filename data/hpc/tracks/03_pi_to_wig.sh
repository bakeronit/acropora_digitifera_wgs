#!/bin/bash

for f in *.pestPG;do
	pop=${f%.pestPG}
	cat $f | awk -v stat="P" -f theta2wig.awk > ${pop}_ThetaD.wig

	wigToBigWig ${pop}_ThetaD.wig adig2_chromsizes.txt ${pop}_ThetaD.bw
done

