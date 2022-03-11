# order A, C, G, T. For example, the first line in the example data file is:
# 20,0,0,0        0,0,0,1 0,0,0,1 0,0,0,1


#%CHROM\t%POS0\t%POS0\t%REF\t%ALT{0}\t%AC\t%AN
# Typical input lines
#
#BLFC01000006.1	76263	76263	G	T	2	150	.	-1	-1	.
#BLFC01000006.1	76266	76266	T	A	149	150	BLFC01000006.1	76266	76266	T,A,A

function printacs(a1,n1,a2,n2)
{
	acs["A"]=0
	acs["C"]=0
	acs["G"]=0
	acs["T"]=0

	for (allele in acs){
		if ( toupper(a1) == allele){
			acs[allele]=n1
		}
		if ( toupper(a2) == allele){
			acs[allele]=n2
		}
	}

	return sprintf("%s,%s,%s,%s",acs["A"],acs["C"],acs["G"],acs["T"])
}

$8=="." {
	print printacs($4,($7-$6),$5,$6), printacs($4,1), printacs($4,1)
}

$8!="." {
	split($11,a,",")
	print printacs($4,($7-$6),$5,$6), printacs(a[2],1), printacs(a[3],1)
}

