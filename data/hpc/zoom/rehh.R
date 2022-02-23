library(tidyverse)
library(rehh)
library(phytools)


position <- 282923
selected_allele <- 1
vcf_file <- "BLFC01000154.1_250001.vcf"
sweep_name <- "BLFC01000154_PXDN"

cat("Reading ",vcf_file," \n")

hap_names2pops <- function(hn){
  pops <- rep("IN",length(hn))
  sos <- grepl("RS[123]",hn)
  nos <- grepl("^AR",hn)
  ins <- (!sos & !nos)
  pops[sos]="SO"
  pops[nos]="NO"
  pops
}

hh <- data2haplohh(hap_file = vcf_file,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

marker_index <- which(positions(hh)==position)

# ---- Furcation
furcation <- calc_furcation(hh,mrk = marker_index)

hn <- hap.names(hh)

newick_r <- as.newick(furcation,
                    allele = selected_allele,
                    side = "right",
                    hap.names = hap.names(hh))

newick_l <- as.newick(furcation,
                      allele = selected_allele,
                      side = "left",
                      hap.names = hap.names(hh))



# ---- Haplotypes
if(marker_index<201){
  s <- 1
  mrk <- marker_index
} else {
  s <- marker_index-200
  mrk <- 201
}
e <- min(marker_index+200,length(positions(hh))+1)

hh_subset <- subset(hh, select.mrk = s:e)






# ---- Trees
hh_gg <- hh_subset@haplo %>% 
  as.data.frame() %>% 
  rownames_to_column("haplotype") %>% 
  pivot_longer(-haplotype,names_to = "position", values_to = "allele") %>% 
  extract(position,into = "pos",regex = "V([0-9]+)",convert = TRUE) %>% 
  filter(allele==1)


library(ape)
tree_l <- ape::read.tree(text = newick_l)
tree_r <- ape::read.tree(text = newick_r)

library(ggtree)

tree_data <- data.frame(tiplab=hn) %>% 
  mutate(pop=hap_names2pops(hn), marker_allele = as.factor(hh_subset@haplo[,mrk]))

tree_data

gtr <- ggtree(tree_r,ladderize = FALSE) %<+% tree_data + 
  geom_tiplab(aes(color=pop),size=2) + 
  geom_tippoint(aes(shape=marker_allele)) +
  geom_facet(panel="Haplotypes",data = hh_gg, geom= geom_point,mapping = aes(x=pos),color="blue",shape="|",size=1)

gtl <- ggtree(tree_l,ladderize = FALSE) %<+% tree_data + 
  geom_tiplab(aes(color=pop),size=2) + 
  geom_tippoint(aes(shape=marker_allele)) +  
  geom_facet(panel="Haplotypes",data = hh_gg, geom= geom_point,mapping = aes(x=pos),color="blue",shape="|",size=1)




njtree <- hh_subset@haplo %>% dist.gene() %>% nj %>% midpoint.root()

source("color_scheme.R")

ggtree(njtree,ladderize = FALSE) %<+% tree_data + 
#  geom_tiplab(aes(label=pop,color=pop),align = TRUE,size=2) +
  geom_tippoint(aes(shape=marker_allele,color=pop),size=1) + 
  geom_facet(panel="Haplotypes",data = hh_gg %>% filter(allele==1), 
             geom= geom_point,mapping = aes(x=pos,color=pop),shape="|",size=1) +
  scale_color_manual(values=myCol) + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    legend.position = "none"
  )

ggsave("tree_haps.png",height = 3,width = 5)

ggtree(njtree,ladderize = FALSE) %<+% tree_data + 
  geom_tiplab(aes(color=pop),align = TRUE,size=1) +
  geom_tippoint(aes(shape=marker_allele,color=pop),size=2) + 
#  geom_facet(panel="Haplotypes",data = hh_gg %>% filter(allele==1), 
#             geom= geom_point,mapping = aes(x=pos,color=pop),shape="|",size=1) +
  scale_color_manual(values=myCol) 

ggsave("tree.png",height = 6,width = 8)

#plot(furcation,hap.names = hap_names2pops(hap.names(hh)),cex=0.5)

