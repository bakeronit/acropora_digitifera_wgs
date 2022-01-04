
echo "BEGIN TRAITS;
	Dimensions NTRAITS=4;
	Format labels=yes missing=? separator=Comma;
	TraitLabels IN NO SO JP;
	Matrix"

cat AllSamplesMitoConsensus.fasta | bioawk -c fastx '$name~/^AI/{printf("%s 1,0,0,0\n", $name)}' 
cat AllSamplesMitoConsensus.fasta | bioawk -c fastx '$name~/^BR/{printf("%s 1,0,0,0\n", $name)}' 
cat AllSamplesMitoConsensus.fasta | bioawk -c fastx '$name~/^AR/{printf("%s 0,1,0,0\n", $name)}' 
cat AllSamplesMitoConsensus.fasta | bioawk -c fastx '$name~/^RS/{printf("%s 0,0,1,0\n", $name)}'
cat AllSamplesMitoConsensus.fasta | bioawk -c fastx '$name~/^DR/{printf("%s 0,0,0,1\n", $name)}'

echo ";"

echo "END;"