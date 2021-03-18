
cal(){
 pop1=$1
 pop2=$2
 grep -v "nan" ${pop1}_${pop2}_per_snp_fst.weir.fst |\
 awk '{
        sum+=$3; sumsq+=$3^2} 
        END {printf "Average fst is %f, with a sd of %f \n", sum/NR, sqrt((sumsq-sum^2/NR)/NR)}'
}

export -f cal

cal inshore northoffshore
cal inshore southoffshore
cal northoffshore southoffshore
