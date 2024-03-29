library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggsci)

samples <- scan("data/hpc/admixture/sample_order.txt", character())
read_Q <- function(filename) {
  tb <- read_table(filename,col_names = FALSE) %>% 
    add_column(samples=samples) %>% 
    mutate(pop=str_split(samples,"_") %>% map_chr(first)) %>% 
    gather(cluster,proportion,-samples,-pop)
  
  tb$samples <- factor(tb$samples, levels = tb %>% 
                         arrange(pop) %>% pull(samples) %>% unique())
  tb$pop <- factor(tb$pop,levels = c("AI","BR","AR","RS1","RS2","RS3"))
  
  tb$cluster <- factor(tb$cluster, levels=rev(unique(tb$cluster)))
  
  tb
}
pop.labs <- c("AI","BR","AR","RS1","RS2","RS3")
names(pop.labs) <- c("AI","BR","AR","RS1","RS2","RS3")

p1 <- ggplot(read_Q("data/hpc/admixture/Adigitifera_ldpruned.2.Q") %>% filter(samples!="BR_5_121_S125_L004"),aes(x=samples)) + 
  geom_col(aes(group=samples,y=proportion,fill=cluster),show.legend = FALSE) + 
  theme_minimal() + 
  facet_grid(~pop,switch = "x",scales = "free",space = "free") +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),strip.text.x = element_blank())+
  labs(y="Ancestry (K=2)", x="") 

ad3<- read_Q("data/hpc/admixture/Adigitifera_ldpruned.3.Q") %>% filter(samples!="BR_5_121_S125_L004")
ad3$pop <- factor(ad3$pop, levels = c("AR","BR","AI","RS1","RS2","RS3"))
ad3$cluster <- factor(ad3$cluster, levels = c("X2","X1","X3"))
p2<-ggplot(ad3,aes(x=samples)) + 
  geom_col(aes(group=samples,y=proportion,fill=cluster),show.legend = FALSE) + 
  theme_minimal(base_size = 12) +
  facet_grid(~pop,switch = "x",scales = "free",space = "free",
             labeller = labeller(pop=pop.labs)) +
  scale_fill_startrek(labels= c("X1"="South Offshore","X2"="North Offshore", "X3" = "Inshore")) +
  scale_y_continuous(expand = c(0,0), position = "right", breaks = c(0,0.5,1)) +
  theme(axis.text.x = element_blank(),panel.grid = element_blank(),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        strip.text.x= element_blank(),
        axis.text.y  = element_text(angle=90),
        axis.title.y.right = element_text(angle=90),
        axis.text.y.right = element_text(hjust=.6))+ 
  labs(y="Ancestry (K=3)", x="") 

p2

ggsave(p2,filename = "Fig.1B.pdf",width = 4,height = 1.35)


pca <- read.table("data/hpc/pca/Adigitifera_smartpca.evec",
                  comment = "#",sep="",
                  col.names =c("sampleid",paste0("PC",1:20),"location")) %>% 
                  mutate(pop=case_when(location %in% c("AI","BR")~"Inshore", 
                                       location == "AR"~"North Offshore", 
                                       location %in% c("RS1","RS2","RS3")~"South Offshore"))

sdev_eval <- sqrt(scan("data/hpc/pca/Adigitifera_smartpca.eval")[1:20])
sum_eval <- sum(scan("data/hpc/pca/Adigitifera_smartpca.eval"))
pve <- sdev_eval/sum_eval*100

pca_plot <- ggplot(pca %>% filter(sampleid!="BR_5_121:BR_5_121_S125_L004"), aes(-PC1,-PC2,color=pop)) +
  geom_point(size=2,alpha=0.8,shape=16)  +
  geom_point(color="black",size=.2) +
  scale_color_startrek() + 
  labs(x=paste0("PC1 (", signif(pve[1],3),"%)"), y=paste0("PC2 (", signif(pve[2],3),"%)")) +
  theme_test(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = c(0.70,0.2),
        legend.text = element_text(size=8),
        legend.key.size = unit(.4, "cm"),
        legend.background = element_rect(fill=alpha('white', 0)),
        #legend.box = "horizontal",
        #legend.direction = "horizontal",
        panel.grid = element_blank()) 
#ggsave(pca_plot,filename = "Fig.1C.pdf",width = 2.5, height = 2.)

d2s_wa <- read_tsv("data/hpc/symbiont/d2s/wa_matrix.txt", col_names = F) %>% 
  column_to_rownames("X1") %>% as.matrix() 
colnames(d2s_wa) <- rownames(d2s_wa)

mds_wa <- d2s_wa %>% cmdscale() %>%  as.data.frame %>% rownames_to_column("id")%>% 
  mutate(location=str_split( string = id,pattern = "_") %>% 
           map_chr(first), pop=case_when(location %in% c("AI","BR")~"inshore",
                                         location%in% c("RS1","RS2","RS3")~"southoffshore", 
                                         location=="AR"~"northoffshore")) %>% mutate(V1=-V1)

mds_plot<- ggscatter(mds_wa, x="V1", y="V2", 
          #color="groups",
          repel=TRUE, 
          size=.5,
          color="pop",
          #color="location",
          ellipse=TRUE,
          #palette = pal_startrek("uniform")(7),
          palette = pal_startrek("uniform")(7),
          ellipse.type="convex",
          ggtheme = theme_test(base_size = 12,base_line_size = NA),legend="none",
          xlab="Dimension 1", ylab="Dimension 2")

#ggsave(mds_plot,filename = "Fig.1D.pdf",width = 2.5,height = 2.2)

plot_grid(pca_plot, mds_plot, nrow = 1, label_size = 12,align = "h", axis = "b", rel_heights = c(0.5,0.5), labels = c("D","E"))
#ggsave("Fig.1DE.pdf",width = 4.8, height = 2.2)
