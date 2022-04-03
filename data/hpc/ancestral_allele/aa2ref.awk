/^#/{print}

!/^#/{
	info=$8
	split(info,aaa,"=")
	aa=aaa[2]
	ref=$4
	alt=$5
	if (aa!=ref){
		if(alt==aa){ # This is by far the most common case
			$4=alt
			$5=ref
			printf("%s",$1)
			for(i=2;i<=NF;i++){
				if(i>9){
					gsub("0","A",$i)
					gsub("1","0",$i)
					gsub("A","1",$i)				
				}
				printf("\t%s",$i)
			}
		printf("\n")
		} else {
			print $0 # In this situation both ref and alt are different from aa. No point swapping then because we don't know which to call "derived"
		}
	} else {
		print $0
	}
}
