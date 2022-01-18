BEGIN {
	min=mv
}

# IN Private
($6==0) && ($7==0){
	if ( $5>=min){
		print
	}
}

($6==1) && ($7==1){
	if ( $5<=(1-min)){
		print
	}
}


# NO Private
($5==0) && ($7==0){
	if ( $6>=min){
		print
	}
}

($5==1) && ($7==1){
	if ( $6<=(1-min)){
		print
	}
}

# SO Private
($5==0) && ($6==0){
	if ( $7>=min){
		print
	}
}

($5==1) && ($6==1){
	if ( $7<=(1-min)){
		print
	}
}
