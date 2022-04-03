# Expected inputs .. [ACGT]
# BLFC01000154.1	816157	816158	G	T	3	150	.	-1	-1	.	32324 133 0.999987 0.000001 0.000000 0.999927 0.000072 
# BLFC01000154.1	816171	816172	A	T	108	150	BLFC01000154.1	816171	816172	A,T,T	32325 316 0.999843 0.000067 0.000001 0.000000 0.999933 

# Expected outputs
# BLFC01000154.1	816158	G
# BLFC01000154.1	816172	T

{
	aa="."
	threshold=0.7 #
	if ( $15 > threshold){
		aa="A"
	} else if ( $16 > threshold){
		aa="C"
	} else if ( $17 > threshold){
		aa="G"
	} else if ( $18 > threshold){
		aa="T"
	}
	printf("%s\t%s\t%s\n",$1,$3,aa)
}