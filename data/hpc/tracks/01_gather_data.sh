cat ../../genome/GCA_014634065.1_Adig_2.0_genomic.fna | bioawk -c fastx '{print $name,length($seq)}' > adig2_chromsizes.txt

ln -s ../thetastat/inshore.pestPG .
ln -s ../thetastat/north_offshore.pestPG .
ln -s ../thetastat/south_offshore.pestPG .

