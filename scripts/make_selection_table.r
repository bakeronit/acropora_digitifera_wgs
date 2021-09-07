## make the tables
library(tidyverse)

pop="inshore"
target_regions <- read_tsv(paste0("data/hpc/selection/temp/",pop,".multiIntersect.filtered.bed")) %>% 
  select(chrom, start,end,num,list)


add_values <- function(chrom, start, end) {
  ihs_file <- paste0("data/hpc/selection/temp/ihs/",chrom,".ihs.out.50bins.norm")
  xpehh_no_file <- paste0("data/hpc/selection/temp/xpehh_no/", chrom, ".xpehh.out.norm")
  xpehh_so_file <- paste0("data/hpc/selection/temp/xpehh_so/", chrom, ".xpehh.out.norm")
  xpnsl_no_file <- paste0("data/hpc/selection/temp/xpnsl_no/", chrom, ".xpnsl.out.norm")
  xpnsl_so_file <- paste0("data/hpc/selection/temp/xpnsl_so/", chrom, ".xpnsl.out.norm")
  
  if(!file.exists(ihs_file)){
    system(paste0("scp zodiac:/home/jc502059/bioprojects/adigitifera_wgs_v2/08.selscan/ihs_out/inshore/",chrom,".ihs.out.50bins.norm ",ihs_file))
  }
  ihs_value <- read_tsv(ihs_file, 
                         col_names = c('chrom','pos','maf','ihha','ihhd','ihs','norm_ihs','sd')) %>% 
    na.omit() %>% filter(pos>=start & pos <end) %>% 
    summarise(prop=length(norm_ihs[abs(norm_ihs)>2])/n(), max_ihs=max(abs(norm_ihs))) %>% 
    add_column(chrom=chrom,start=start,end=end)
  
  if(!file.exists(xpehh_no_file)){
    system(paste0("scp zodiac:/home/jc502059/bioprojects/adigitifera_wgs_v2/08.selscan/xpehh_out/inshore_vs_northoffshore/", chrom, ".xpehh.out.norm ",xpehh_no_file))
  }
  xpehh_no_value <- read_tsv(xpehh_no_file) %>% na.omit %>% filter(pos>=start & pos<end) %>% 
    summarise(prop_xpehh_no=length(crit[crit==1])/n(), max_xpehh_no = max(normxpehh))
  
  
  if(!file.exists(xpehh_so_file)){
    system(paste0("scp zodiac:/home/jc502059/bioprojects/adigitifera_wgs_v2/08.selscan/xpehh_out/inshore_vs_southoffshore/", chrom, ".xpehh.out.norm ",xpehh_so_file))
  }
  xpehh_so_value <- read_tsv(xpehh_so_file) %>% na.omit %>% filter(pos>=start & pos<end) %>% 
    summarise(prop_xpehh_so=length(crit[crit==1])/n(), max_xpehh_so = max(normxpehh))
  
  if(!file.exists(xpnsl_no_file)){
    system(paste0("scp zodiac:/home/jc502059/bioprojects/adigitifera_wgs_v2/08.selscan/xpnsl_out/inshore_vs_northoffshore/", chrom, ".xpnsl.out.norm ",xpnsl_no_file))
  }
  xpnsl_no_value <- read_tsv(xpnsl_no_file) %>% na.omit %>% filter(pos>=start & pos<end) %>% 
    summarise(prop_xpnsl_no=length(crit[crit==1])/n(), max_xpnsl_no = max(normxpehh))
  
  if(!file.exists(xpnsl_so_file)){
    system(paste0("scp zodiac:/home/jc502059/bioprojects/adigitifera_wgs_v2/08.selscan/xpnsl_out/inshore_vs_southoffshore/", chrom, ".xpnsl.out.norm ",xpnsl_so_file))
  }
  xpnsl_so_value <- read_tsv(xpnsl_so_file) %>% na.omit %>% filter(pos>=start & pos<end) %>% 
    summarise(prop_xpnsl_so=length(crit[crit==1])/n(), max_xpnsl_so = max(normxpehh))
  
  cbind(ihs_value,xpehh_no_value, xpnsl_no_value, xpehh_so_value, xpnsl_so_value)
}


add_values("BLFC01000008.1",1150001,1200001)

tb <- pmap_df(list(target_regions$chrom, target_regions$start, target_regions$end), ~add_values(..1,..2,..3)) %>% as_tibble


bl6_cols <- c("qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
blastp <- read_tsv("data/hpc/annotation/Adigitifera_uniprot_blastp.outfmt6", col_names = bl6_cols) %>% 
  select(adi_id=qaccver, uniprot_id=saccver,evalue) %>% add_column(method="blastp") %>% mutate(geneid=str_remove(adi_id,'\\.t\\d'))

blastx <- read_tsv("data/hpc/annotation/Adigitifera_uniprot_blastx.outfmt6", col_names = bl6_cols) %>% 
  select(adi_id=qaccver, uniprot_id=saccver,evalue) %>% add_column(method="blastx") %>% mutate(geneid=str_remove(adi_id,'\\.t\\d'))

uniprotkb_tab <- read_tsv("data/hpc/annotation/uniprot-yourlist.tab") %>% select(-Status)

uniprot_gene_annot <- rbind(blastp,blastx) %>% group_by(geneid) %>% dplyr::top_n(1,dplyr::desc(evalue)) %>% group_by(geneid,evalue) %>% filter(row_number()==1) %>% 
  left_join(uniprotkb_tab,by=c("uniprot_id"="Entry")) %>%
  ungroup %>% 
  select(geneid,uniprot_id,entryname="Entry name",
         genename="Gene names",go="Gene ontology IDs",kegg="Cross-reference (KEGG)",protein="Protein names") %>% 
  mutate(genename=ifelse(is.na(genename),entryname,genename))


genes <- read_tsv("data/hpc/selection/temp/inshore.test.genes.bed", col_names = c("chrom","start","end","geneid"))
snp_count <- read_tsv("data/hpc/selection/temp/snp_counts.txt", col_names = c("chrom","start","end","n_snp"))
final_tb <- target_regions  %>% mutate(chrom=str_remove(chrom,"\\.1")) %>% 
  left_join(genes) %>% left_join(uniprot_gene_annot) %>% 
  mutate(entryname=ifelse(is.na(entryname),geneid,entryname)) %>% 
  group_by(chrom, start, end, num, list) %>% 
  summarise(genes = paste(entryname, collapse = "; ")) %>% ungroup %>% 
  left_join(tb %>% mutate(chrom=str_remove(chrom,"\\.1"))) %>%
  mutate(genes=ifelse(genes=="NA","no genes in this region", genes)) %>% 
  arrange(desc(num),desc(prop)) %>% left_join(snp_count %>% mutate(chrom=str_remove(chrom,"\\.1")))

write_tsv(final_tb, "selection_table.txt")

