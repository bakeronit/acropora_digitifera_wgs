# AA Code [ACDEFGHIKLMNPQRSTVWY]

/^>/ {
	print
}

!/^>/ {
	gsub(/[^ACDEFGHIKLMNPQRSTVWY]/,"")
	print
}