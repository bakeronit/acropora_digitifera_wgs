# The first three required BED fields are:

# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
# The 9 additional optional BED fields are:

# name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
# score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
# shade	 	 	 	 	 	 	 	 	 
# score in range  	≤ 166	167-277	278-388	389-499	500-611	612-722	723-833	834-944	≥ 945
# strand - Defines the strand. Either "." (=no strand) or "+" or "-".
# thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
# thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
# itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
# blockCount - The number of blocks (exons) in the BED line.
# blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
# blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.


#BLFC01000154.1	AUGUSTUS	gene	223522	252915	0.11	-	.	ID=adig_s0150.g21
#BLFC01000154.1	AUGUSTUS	gene	261498	262495	0.17	-	.	ID=adig_s0150.g22
#BLFC01000154.1	AUGUSTUS	gene	262827	267715	0.03	-	.	ID=adig_s0150.g23
#BLFC01000154.1	AUGUSTUS	gene	275513	284565	0.17	-	.	ID=adig_s0150.g24
#BLFC01000154.1	AUGUSTUS	gene	296962	300915	0.05	-	.	ID=adig_s0150.g25


($3=="gene") && ($4>150000) && ($5<400000){
	OFS="\t"

	match($9,"adig_s0[^\;]*")

	name=substr($9,RSTART,RLENGTH)

	is_highlighted="no"


	display_name=""
	display_color="169,169,169"

	if ( name=="adig_s0150.g18"  ){
		display_name="PXDN s0150.g18"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g19"  ){
		display_name="PXDN s0150.g19"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g20"  ){
		display_name="PXDN s0150.g20"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g21"  ){
		display_name="PXDN s0150.g21"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g22"  ){
		display_name="s0150.g22"
		is_highlighted="no"
	}

	if ( name=="adig_s0150.g23"  ){
		display_name="PERC s0150.g23"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g24"  ){
		display_name="PXDN s0150.g24"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g25"  ){
		display_name="PXDN s0150.g25"
		display_color="0,150,255"
		is_highlighted="yes"
	}

	if ( name=="adig_s0150.g26"  ){
		display_name="PXDN s0150.g26"
		display_color="0,150,255"
		is_highlighted="yes"
	}






	if ( is_highlighted=="yes"){
		print $1,$4,$5,display_name,0,$7,$4,$5,display_color > "px.bed"
	} else {
		print $1,$4,$5,display_name,0,$7,$4,$5,display_color > "nopx.bed"
	}

	print $1,$4,$5,display_name,0,$7,$4,$5,display_color

}