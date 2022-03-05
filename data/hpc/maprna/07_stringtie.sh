for f in *_BLFC01000154.1.bam;do
	sample=${f%_BLFC01000154.1.bam}
	stringtie $f -G ../../genome/adig-v2-ncbi.gff --conservative -p 60 -A ${sample}_geneab.out -o ${sample}.gtf
done

ls *.gtf > mergelist.txt

stringtie --merge -p 40 -G ../../genome/adig-v2-ncbi.gff -o BLFC01000154.1_stringtie_merged.gtf mergelist.txt
