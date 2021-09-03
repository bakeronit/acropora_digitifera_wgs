bl6_cols <- c("qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
blastp <- read_tsv("data/hpc/annotation/Adigitifera_uniprot_blastp.outfmt6", col_names = bl6_cols) %>% 
  dplyr::select(adi_id=qaccver, uniprot_id=saccver,evalue) %>% add_column(method="blastp") %>% mutate(geneid=str_remove(adi_id,'\\.t\\d'))

blastx <- read_tsv("data/hpc/annotation/Adigitifera_uniprot_blastx.outfmt6", col_names = bl6_cols) %>% 
  dplyr::select(adi_id=qaccver, uniprot_id=saccver,evalue) %>% add_column(method="blastx") %>% mutate(geneid=str_remove(adi_id,'\\.t\\d'))

uniprotkb_tab <- read_tsv("data/hpc/annotation/uniprot-yourlist.tab") %>% dplyr::select(-Status)

uniprot_gene_annot <- rbind(blastp,blastx) %>% group_by(geneid) %>% dplyr::top_n(1,dplyr::desc(evalue)) %>% group_by(geneid,evalue) %>% filter(row_number()==1) %>% 
  left_join(uniprotkb_tab,by=c("uniprot_id"="Entry")) %>%
  ungroup %>% 
  dplyr::select(geneid,uniprot_id,entryname="Entry name",
         genename="Gene names",go="Gene ontology IDs",kegg="Cross-reference (KEGG)") %>% mutate(genename=ifelse(is.na(genename),entryname,genename))