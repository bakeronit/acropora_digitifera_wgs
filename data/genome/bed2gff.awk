
# The first three fields in each feature line are required:

# 1 chrom - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be given with or without the 'chr' prefix.
# 2 chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
# 3 chromEnd - End position of the feature in standard chromosomal coordinates
# 4 name - Label to be displayed under the feature, if turned on in "Configure this page".
# 5 score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
# 6 strand - defined as + (forward) or - (reverse).
# thickStart - coordinate at which to start drawing the feature as a solid rectangle
# thickEnd - coordinate at which to stop drawing the feature as a solid rectangle
# itemRgb - an RGB colour value (e.g. 0,0,255). Only used if there is a track line with the value of itemRgb set to "on" (case-insensitive).
# blockCount - the number of sub-elements (e.g. exons) within the feature
# blockSizes - the size of these sub-elements
# blockStarts - the start coordinate of each sub-element

# sc0000260_arrow_pilon	AUGUSTUS	gene	18066	21459	0.36	+	.	ID=adig_s0260.g1

#BLFC01000061.1	18065	21459	transcript	0	+	ID=adig_s0260.g1.t1;Parent=adig_s0260.g1

BEGIN {
	OFS="\t"
}

{
	print $1,"AUGUSTUS",$4,$2+1,$3,0,$6,$5,$7
}