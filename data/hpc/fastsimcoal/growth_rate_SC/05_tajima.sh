cd 3.out.growth_rate_SC

mkdir -p tajima

for f in *.vcf.gz;do
    echo $f
    for pop in IN NO SO;do
        echo $pop
        bcftools view -S ${pop}.pop ${f} |\
         vk tajima 10000 2000 - |\
         sed '1d' |\
         awk '{print $1"\t"$2+1"\t"$3"\t"$6}' > tajima/${f}_${pop}_1based.td
    done
done
