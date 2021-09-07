# this script is to make cms like window plot for any window we provided.
library(tidyverse)
library(cowplot)
library(ggrepel)
library(gridExtra)

ihs <- read_tsv("data/hpc/selection/windows/multiIntersect.ihs.norm", 
                col_names = c("scaff","pos","frq","ihh1","ihh2","unstand_ihs","value","sd")) %>% 
  select(scaff, pos, value) %>% add_column(type="ihs")

xpehh <- read_tsv("data/hpc/selection/windows/multiIntersect.xpehh.norm", 
                  col_names = c("scaff","pos","gpos","p1","ihha","ihhb","p2","uxpehh","value","sd")) %>% 
  select(scaff,pos,value) %>% add_column(type="xpehh")

xpnsl <- read_tsv("data/hpc/selection/windows/multiIntersect.xpnsl.norm", 
                  col_names = c("scaff","pos","gpos","p1","ihha","ihhb","p2","uxpnsl","value","sd")) %>% 
  select(scaff,pos,value) %>% add_column(type="xpnsl")

fst <- read_tsv("data/hpc/selection/windows/multiIntersect.weir.fst",
                col_names = c("scaff","pos","value")) %>% 
  add_column(type="fst")

dihh <- read_tsv("data/hpc/selection/windows/multiIntersect.deltaiHH.txt",
                 col_names = c("scaff","pos","value")) %>% add_column(type="deltaiHH")
ddaf <- read_tsv("data/hpc/selection/windows/multiIntersect.deltaDAF.txt",
                 col_names = c("scaff","pos","value")) %>% add_column(type="deltaDAF")

bed <- read_tsv("data/hpc/selection/results/inshore.multiIntersect.bed") %>% 
  select(chrom,start,end) %>% rename("scaff"=chrom)

allgenes <- read_tsv("data/hpc/selection/results/allgenes.bed", 
                     col_names = c("scaff","gstart","gend","adi_id"))

window_genes <-  bed %>% left_join(allgenes %>% mutate(scaff=paste0(scaff,".1"))) %>% 
  filter(gstart >= start & gend <= end) %>% 
  mutate(location=paste(scaff,paste(start,end,sep="-"),sep = ":")) %>% 
  select(location, gstart,gend, adi_id)

uniprot_gene_annot <- read_rds("cache/uniprot_gene_annot.rds")
window_genes <- window_genes %>% left_join(uniprot_gene_annot,by=c("adi_id"="geneid")) %>% 
  select(location, gstart, gend, adi_id, entryname) %>% 
  mutate(entryname=ifelse(is.na(entryname),adi_id,entryname))

window_data <- rbind(ihs, xpehh, xpnsl, fst, dihh, ddaf) %>% left_join(bed) %>%  
  filter(pos>=start & pos<=end) %>% 
  mutate(location=paste(scaff,paste(start,end,sep="-"),sep = ":")) %>% 
  select(location, pos, value, type)




plot_one_window <- function(region, scale=1e+6,unit="Mb") {
  scaff <- str_extract(region,"BLFC\\d+\\.\\d")
  begin <- str_match(region,"\\:(\\d+)\\-")[,2] %>% as.integer
  end <-str_match(region,"\\-(\\d+)")[,2] %>% as.integer
  
  p_ihs <- window_data %>% filter(location==region & type=="ihs") %>% 
    ggplot(aes(x=pos/scale,y=value)) + geom_point(size=.1,color="darkblue") + theme_bw() + 
    theme(panel.grid = element_blank(),axis.title.x = element_blank()) + labs(y="iHS") + xlim(begin/scale,end/scale) 
  
  p_xpehh <- window_data %>% filter(location==region & type=="xpehh") %>% 
    ggplot(aes(x=pos/scale,y=value)) + geom_point(size=.1,color="darkblue") + theme_bw() + 
    theme(panel.grid = element_blank(),axis.title.x = element_blank()) + labs(y="XP-EHH") + xlim(begin/scale,end/scale)  
  
  p_xpnsl <- window_data %>% filter(location==region & type=="xpnsl") %>% 
    ggplot(aes(x=pos/scale,y=value)) + geom_point(size=.1,color="darkblue") + theme_bw() + 
    theme(panel.grid = element_blank(),axis.title.x = element_blank()) + labs(y="XP-nSL") + xlim(begin/scale,end/scale) 
  
  p_dihh <- window_data %>% filter(location==region & type=="deltaiHH") %>% 
    ggplot(aes(x=pos/scale,y=value)) + geom_point(size=.1,color="darkblue") + theme_bw() + 
    theme(panel.grid = element_blank(),axis.title.x = element_blank()) + labs(y=expression(paste(Delta,"iHH"))) + xlim(begin/scale,end/scale) 
  
  p_ddaf <- window_data %>% filter(location==region & type=="deltaDAF") %>% 
    ggplot(aes(x=pos/scale,y=value)) + geom_point(size=.1,color="darkblue") + theme_bw() + 
    theme(panel.grid = element_blank(),axis.title.x = element_blank()) + labs(y=expression(paste(Delta,"DAF"))) + xlim(begin/scale,end/scale)  
  
  p_fst <- window_data %>% filter(location==region & type=="fst") %>% 
    ggplot(aes(x=pos/scale,y=value)) + geom_point(size=.1,color="darkblue") + theme_bw() + 
    theme(panel.grid = element_blank()) + labs(x=paste("Position on scaffold ",scaff, "in",unit),y="Fst") + xlim(begin/scale,end/scale) 
  
  p_genes <- window_genes %>% filter(location==region) %>% 
    ggplot() + geom_segment(aes(x=gstart/scale, xend=gend/scale, y=1,yend=1),size=2,linetype=1) + 
    geom_label(aes(x=(gstart/scale+gend/scale)/2,label=entryname,y=1),size=2) + ylim(0,1.2)  + xlim(begin/scale,end/scale) +
    theme_classic() + theme(axis.line = element_blank(),
                            axis.title = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text = element_blank()) 

  plot_grid(p_ihs,p_xpehh,p_xpnsl,p_dihh,p_ddaf,p_fst,p_genes,ncol = 1)
}


plot_one_window("BLFC01000770.1:1200001-1250001")

#TODO: how to plot multiple plot in one pdf???
pdf("test.pdf", onefile = TRUE)
for (i in 1:9){
 do.call("grid.arrange", plot_one_window(window_data$location %>% unique %>% nth(i)))
}
dev.off()





                          