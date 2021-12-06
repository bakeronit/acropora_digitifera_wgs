BEGIN {
	current_chr=""
}

NR>1 {
	if ( current_chr!=$2 ){
		current_chr=$2
		printf("variableStep chrom=%s\n",current_chr)
	} else {
		if ( stat == "P"){
			printf("%i\t%s\n",$3,$5)		
		} else {
			printf("%i\t%s\n",$3,$9)
		}
	}
}
