{
	fst1=$5
	fst2=$10
	fst3=$15
	tinno=-log(1-fst1)
	tinso=-log(1-fst2)
	tnoso=-log(1-fst3)

	pbsin=(tinso+tinso-tnoso)/2
	pbsso=(tnoso+tinso-tinno)/2
	pbsno=(tnoso+tinno-tinso)/2

	print $1,$2,pbsin,pbsno,pbsso

}