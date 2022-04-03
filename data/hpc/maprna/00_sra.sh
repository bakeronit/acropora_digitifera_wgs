while read acc;do
  echo "Downloading $acc"
  singularity run $SING/sra-tools-2.11.0.sif prefetch ${acc}
  singularity run $SING/sra-tools-2.11.0.sif fasterq-dump --split-files ${acc}/${acc}.sra
  echo "Zipping $acc"
  gzip ${acc}/${acc}_1.fastq
  gzip ${acc}/${acc}_2.fastq
done < <(cat SRR_Acc_List*.txt)
