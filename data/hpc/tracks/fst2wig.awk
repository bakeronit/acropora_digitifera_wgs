BEGIN {
	current_chr=""
	col=pop+8
}

NR>1 {
	if ( current_chr!=$2 ){
		current_chr=$2
		printf("variableStep chrom=%s\n",current_chr)
	} else {
		printf("%i\t%s\n",$3,$col)
	}
}
