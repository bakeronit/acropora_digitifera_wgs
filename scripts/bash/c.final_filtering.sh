#1. filtering by depth, genotype missingness etc.
vcftools --gzvcf Adigi.indel5bp_snp_hardfilter_passed_biallelic_mdust.vcf.gz \
    --max-missing 0.9 --minQ 30 \
    --min-meanDP 10 --max-meanDP 32 \
    --minDP 3 --minGQ 20 --remove-filtered-geno-all \
    --recode --recode-INFO-all \
    --stdout | bgzip > Adigi.v2.DPg90gdp3gq30.vcf.gz

bcftools stats --threads 20 Adigi.v2.DPg90gdp3gq30.vcf.gz > 6.Adigi.v2.DPg90gdp3gq30.stats 


#2. filtering snps with Inbreeding Coefficient < -0.05
gatk VariantFiltration \
    -V  Adigi.v2.DPg90gdp3gq30.vcf.gz\
    -filter "InbreedingCoeff < -0.05" --filter-name "ICF0.05" \
    -O Adigi.v2.DPg90gdp3gq30.Fis0.05.vcf.gz

gatk SelectVariants \
    -V Adigi.v2.DPg90gdp3gq30.Fis0.05.vcf.gz \
    -O Adigi.v2.DPg90gdp3gq30.Fis0.05_pass.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    --exclude-filtered true

bcftools stats --threads 20 Adigi.v2.DPg90gdp3gq30.Fis0.05_pass.vcf.gz > 7.Adigi.v2.DPg90gdp3gq30_Fis0.05_pass.stats

# remove monomorphic sites
bcftools filter -e 'AC=0 || AC==AN' --threads 20 Adigi.v2.DPg90gdp3gq30.Fis0.05_pass.vcf.gz |bgzip > Adigi.v2.filtered.vcf.gz
bcftools stats --threads 20 Adigi.v2.filtered.vcf.gz > Adigi.v2.filtered.stats
