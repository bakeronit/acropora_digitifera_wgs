library(rehh)
library(vcfR)
library(cowplot)

locus <- read_tsv("data/hpc/selection/vcf_files/locus.txt", col_names = "id")

make_ehh_plot <- function(loci) {

  scaff <- str_split(loci,":")[[1]][1]
  pos <- str_split(loci,":")[[1]][2] %>% as.numeric()
  print(scaff)
  hh <- data2haplohh(hap_file = paste0("data/hpc/selection/vcf_files/", scaff, "_withid.vcf"),
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
  hh_offshore <- data2haplohh(hap_file = paste0("data/hpc/selection/vcf_files/", scaff, "_offshore_withid.vcf"),
                              polarize_vcf = FALSE,
                              vcf_reader = "vcfR")
  
  res <- calc_ehh(hh, mrk = loci, include_nhaplo = FALSE)
  res_offshore <- calc_ehh(hh_offshore, mrk = loci, include_nhaplo = FALSE)
  
  pehh <- res$ehh %>% as_tibble() %>% pivot_longer(-POSITION) %>% 
    ggplot() + geom_line(aes(x=(POSITION-pos)/1e+6,y=value, color=name),size=1.5) +
    theme_classic() + ggsci::scale_color_aaas() +
    theme(legend.position = "none") + 
    labs(x="Position(Mb)",y="Extended Haplotype Homozygosity",title = paste("EHH around", loci,"in inshore"))
  
  pehh2 <- res_offshore$ehh %>% as_tibble() %>% pivot_longer(-POSITION) %>% 
    ggplot() + geom_line(aes(x=(POSITION-pos)/1e+6,y=value, color=name),size=1.5) +
    theme_classic() + ggsci::scale_color_aaas(name="",label=c("Ancestral","Derived")) +
    labs(x="Position(Mb)",y="Extended Haplotype Homozygosity",title = paste("EHH around", loci,"in offshore"))
  
  plot_grid(pehh, pehh2)
}

make_ehh_plot(locus$id[1])



hh <- data2haplohh(hap_file = paste0("data/hpc/selection/vcf_files/BLFC01000156_offshore_withid.vcf"),
                   polarize_vcf = FALSE,
                   vcf_reader = "vcfR")
res <- calc_ehh(hh, mrk = "BLFC01000156.1:905489" , include_nhaplo = TRUE)
plot(res)
