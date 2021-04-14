# The goal here is to simulate the MyBaits V4 target-enrichment process.
# In addition we use this step to do some filtering to ensure we are working 
# With high quality unique probe matches only
# 
# Since DNA is sheared to 400-500bp the assembled contigs from 
# target enrichment tend to be around 800-1000bp long .. ie roughly 500bp each side of the probe. 
# To simulate this we simply take the central point of the alignment 
# And then create an interal +/- 500bp around it. 
#
#
#
#
# samtools view -F 2308 means
# - Exclude;
#   Unmapped reads
#	Non-primary alignments
#	Supplementary alignments
#

for f in *_probes_clip.sam;do
	echo $f
	cat $f | samtools view -F 2308 | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$3,$4-400,$4+600,$1,$5)}' | awk '$2>0 && $3>0' \
	> ${f/.sam/.bed}
done



