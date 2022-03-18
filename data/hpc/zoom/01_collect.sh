cat ../../genome/adig-v2-ncbi.gff | grep 'BLFC01000154.1' | awk '$3=="gene"{print}' > BLFC01000154.gff

cat BLFC01000154.gff | awk -f s0150.gff2bed.awk

cat ../../genome/adig-v2-ncbi.gff | grep 'BLFC01000600.1' | awk '$3=="gene"{print}' > BLFC01000600.gff

cat BLFC01000600.gff | awk -f s0005.gff2bed.awk
