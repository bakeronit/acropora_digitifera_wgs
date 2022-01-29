cat ../../genome/adig-v2-ncbinames.gff | grep 'BLFC01000154.1' | awk '$3=="gene"{print}' > BLFC01000154.gff

cat BLFC01000154.gff | awk -f gff2bed.awk