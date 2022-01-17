#module load kraken
#module load parallel

make_kraken_report(){
	db=$1
  f=$(basename $2);

  r1=$2;
  r2=${r1/R1/R2}

  echo $r1 $r2

  kraken_raw=genome_kraken_out/${f%_R1.fastq}.kraken

  report_out=genome_kraken_reports/${f%R1.fastq}krakenrep.txt
  mpa_out=genome_kraken_mpa/${f%_R1.fastq}.mpa
  printf "Processing %s and %s into %s\n" $r1 $r2 $kraken_raw
  kraken --db $db --threads 8 --fastq-input --paired $r1 $r2 > ${kraken_raw} 

  printf "Running kraken report\n"
  kraken-report --db $db ${kraken_raw} > ${report_out}
  printf "Running kraken-mpa-report\n"
  kraken-mpa-report --db $db ${kraken_raw} > ${mpa_out}
}

export -f make_kraken_report

parallel -j 24 make_kraken_report /fast/shared/kraken_coral_symbionts2 {} ::: $(ls reads/*R1.fastq)
