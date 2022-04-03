# Instructions obtained from https://github.com/ComparativeGenomicsToolkit/hal/blob/chaining-doc/doc/chaining-mapping.md

# Note The UCSC liftOver program considered the target to be the source, so the these instructions create alignments where the target is the source genome.
#
# We are mapping oist (source) onto ncbi (destination)
# So oist is the target and ncbi is the (query)

# First make 2bit files for each of the genome references in our hal file

singularity run -B $(pwd):/data ../hpc/ancestral_allele/cactus_v2.0.5.sif hal2fasta adig.hal oist | faToTwoBit stdin oist.2bit
singularity run -B $(pwd):/data ../hpc/ancestral_allele/cactus_v2.0.5.sif hal2fasta adig.hal ncbi | faToTwoBit stdin ncbi.2bit

# Now export a bed for our source genome sequences
#
singularity run -B $(pwd):/data ../hpc/ancestral_allele/cactus_v2.0.5.sif halStats --bedSequences ncbi /data/adig_amil.hal > ncbi.bed

# Now perform halLiftOver
#
singularity run -B $(pwd):/data ../hpc/ancestral_allele/cactus_v2.0.5.sif halLiftover --outPSL /data/adig_amil.hal ncbi /data/ncbi.bed oist /dev/stdout | pslPosTarget stdin ncbi-oist.psl

# Use UCSC Tools to do the rest. 
# As in this doc http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver
#
axtChain -psl -linearGap=loose ncbi-oist.psl oist.2bit ncbi.2bit ncbi-oist.chain

chainMergeSort ncbi-oist.chain | chainSplit chain stdin

cd chain 
mkdir ../net ../over

for f in *.chain;do
	scaff=${f%.chain}
	chainNet ${scaff}.chain ../oist_lengths.txt ../ncbi_lengths.txt ../net/${scaff}.net /dev/null
	netChainSubset ../net/${scaff}.net ${scaff}.chain ../over/${scaff}.chain
	chainFilter -strand="+" ../over/${scaff}.chain > ../over/${scaff}_pos.chain
done

cd ..

mkdir beds 
while read line;do 
	echo $line; 
	grep $line oist_full.bed > beds/${line}.bed;
done< <(cut -f 1 oist_lengths.txt)

#chainMergeSort over/*.chain > adig_oist2ncbi.chain

cd beds
for bf in *pilon.bed;do
	echo $bf;
	# UCSC liftOver command
	liftOver $bf ../over/${bf%.bed}_pos.chain ${bf%.bed}_converted.bed ${bf%.bed}_unmapped
done

# Consolidate converted and also cleanup, ensuring mapped genes only on good contig pairs
#
rm adig-v2-ncbi.bed
tail -n+2 good_pairs.tsv | awk '{cmd=sprintf("grep %s beds/%s_converted.bed >> adig-v2-ncbi.bed",$2,$1); system(cmd);}'
cat beds/*_unmapped | awk '{print $7}' | grep 'adig' | tr ';' '\n' | awk -F '=' '{print $2}' | sed -E 's/(.*g[0-9]+).t[0-9].*/\1/' | sort -u > unmapped_genes.txt

rm unmapped_gene_patterns.txt
while read gline;do
	g1=$(printf "%s$" $gline)
	echo $g1 >> unmapped_gene_patterns.txt
	g2=$(printf "%s.t" $gline)
	echo $g2 >> unmapped_gene_patterns.txt	
	echo $g1 $g2;
done < unmapped_genes.txt

cat adig-v2-ncbi.bed | sort | uniq > adig-v2-ncbi_su.bed

grep -f unmapped_gene_patterns.txt -v adig-v2-ncbi_su.bed | awk -f bed2gff >  adig-v2-ncbi.gff

