
# This is what the vcf header looks like
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AI_1_001_S102_L004	AI_1_008_merged	AI_1_021_S97_L004	AI_1_022_S98_L004	AI_1_023_S101_L004	AI_1_025_S99_L004	AI_2_036_S105_L004	AI_2_041_S107_L004	AI_2_043_S106_L004	AI_2_136_merged	AI_2_151_S103_L004	AI_3_047_S110_L004	AI_3_060_S109_L004	AI_3_063_S108_L004	AI_3_071_S111_L004	AR_125_374_merged	AR_125_377_S131_L004	AR_125_385_S127_L004	AR_125_388_S128_L004	AR_125_392_merged	AR_128_316_merged	AR_128_318_S135_L004	AR_128_326_merged	AR_128_328_merged	AR_128_336_merged	AR_132_154_merged	AR_132_162_S140_L004	AR_132_170_merged	AR_132_173_S141_L004	AR_132_178_merged	AR_133_341_S143_L004	AR_133_343_merged	AR_133_346_merged	AR_133_354_merged	AR_133_357_S142_L004	BR_4_077_merged	BR_4_078_S116_L004	BR_4_081_S118_L004	BR_4_082_S117_L004	BR_4_087_S112_L004	BR_4_088_merged	BR_4_091_S113_L004	BR_4_100_S119_L004	BR_5_112_S122_L004	BR_5_114_S126_L004	BR_5_121_S125_L004	BR_5_123_S121_L004	BR_5_124_S124_L004	BR_5_129_S120_L004	BR_5_133_S123_L004	RS1_2_417_merged	RS1_2_422_merged	RS1_M11_820_merged	RS1_M11_840_merged	RS1_M12_808_merged	RS1_M12_817_S151_L004	RS1_S_314_S154_L004	RS1_S_321_merged	RS2_2_256_S155_L004	RS2_C11_769_merged	RS2_C11_784_S156_L004	RS2_C13_704_merged	RS2_C13_706_S159_L004	RS2_C13_721_S158_L004	RS2_C20_283_merged	RS2_S_734_S162_L004	RS2_S_737_merged	RS3_1_184_merged	RS3_1_185_S166_L004	RS3_1_191_merged	RS3_1_207_merged	RS3_S_215_merged	RS3_S_232_merged	RS3_S_246_merged
#BLFC01000004.1	8994	.	A	G	.	PASS	.	GT	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0

BEGIN{
print("##fileformat=VCFv4.1")
print("##FILTER=<ID=PASS,Description=\"All filters passed\">")
print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">")
}


/Chrom/{
	printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for (i = 5; i <= NF; i+=2){
		if($i~/A_1/){
			pop="IN"
		}
		if($i~/A_2/){
			pop="NO"
		}
		if($i~/A_3/){
			pop="SO"
		}

		printf("\t%s_%s",pop,(i+1)/2)
	}
	printf("\n")
}

$1 != "Chrom" {
	if (NF!=0){
		printf("contig%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT",$1,$2,$3,$4)
		for (i = 5; i <= NF; i+=2){
        	printf("\t%s|%s",$i,$(i+1))
		}
    	printf("\n")
    }
}