## to plot figure 2

### smc plot
mu <- 1.2e-8 
gen <- 5
#dtime = data.frame(name=c("LGM","Holocene"), start=c(1.2e+4,0), end=c(1.1e+5,1.2e+4))

xbreaks <- c(1e+2*1:9, 1e+3*1:9 , 1e+4*1:9 , 1e+5*1:9 , 1e+6*1:9 , 1e+7)
xlabels <- as.character(xbreaks)
xlabels[!(xbreaks%in%c(1e+3,1e+4,1e+5,1e+6,1e+7))] <- ''
xlabels[xbreaks%in%c(1e+3,1e+4,1e+5,1e+6,1e+7)] <-c("1kya","10kya","100kya","1mya","10mya")

ybreaks <- c(1e+3*2:9 , 1e+4*1:9 , 1e+5*1:9 , 1e+6*1:2)
ylabels <- as.character(ybreaks)
ylabels[!(ybreaks%in%c(1e+4,1e+5,1e+6))] <- ''
ylabels[ybreaks%in%c(1e+4,1e+5,1e+6)] <- c("1","10","100")
#smc <- read_csv("data/hpc/demography/estimate_em50_cubic_k10/em50_cubic_k10.csv")
smc_main <- read_csv("data/hpc/demography/estimate_em50_k40/em50_k40.csv") %>% select(label,x,y) %>% add_column(type="main")

read_smc_range <- function(filename) {
  mu<- filename %>% basename %>% gsub(pattern = "_main_g[3,5,7].csv",replacement = "")
  g<-filename %>% basename %>% gsub(pattern = "mu[1,2,3]_main_(g[3,5,7]).csv",replacement = "\\1")
  read_csv(filename) %>% add_column(mu=mu,g=g)
}

all_smc <- list.files(path = "data/hpc/demography/smc++_new_results", 
                      pattern="mu[1,2,3]_main_g[3,5,7].csv",recursive = T,
                      full.names = T) %>% map_df(read_smc_range) %>% group_by(label,mu,g) %>%  
  slice_min(n=1,y) %>% slice_min(n=1,x) %>% ungroup %>% group_by(label) %>% 
  dplyr::summarise(min_x=min(x),max_x=max(x))

split_in_no <- read_csv("data/hpc/demography/split_em50_cubic_k10/inshore_northoffshore.split.csv")
split_in_so <- read_csv("data/hpc/demography/split_em50_cubic_k10/inshore_southoffshore.split.csv")
split_no_so <- read_csv("data/hpc/demography/split_em50_cubic_k10/northoffshore_southoffshore.split.csv")

split_df <- rbind(split_in_no %>% group_by(label) %>% dplyr::summarise(max=max(x)) %>% add_column(type="split. inshore vs north offshore") %>% slice_min(n=1,max),
split_in_so %>% group_by(label) %>% dplyr::summarise(max=max(x)) %>% add_column(type="split. inshore vs south offshore") %>% slice_min(n=1,max),
split_no_so %>% group_by(label) %>% dplyr::summarise(max=max(x)) %>% add_column(type="split. north offshore vs south offshore") %>% slice_min(n=1,max))


smc_bootstrap_plot<-ggplot() + geom_line(data=smc_main,aes(x=x,y=y,color=label),size=1) +
  scale_color_startrek(labels = c("Inshore", "North Offshore","South Offshore")) +
  #ggnewscale::new_scale_color()+
  #geom_vline(data=split_df, aes(xintercept=max,color=type),linetype="dashed",size=1) + 
  #scale_color_manual(values=c( "#944A6D","#A86500", "#70A36D")) +
  scale_x_log10(breaks=c(1e+3,1e+4,1e+5,1e+6), labels=c("1kya","10kya","100kya","1mya"), limits=c(1e+3,1.3e+6),expand = c(0, 0))+
  scale_y_log10(breaks=ybreaks,labels=ylabels,limits=c(2e+3,2e+6)) + 
  theme_cowplot() + 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.75,0.9), legend.title = element_blank()) +
  labs(x=expression(paste("Years Ago (g=5, ",mu, "=1.20e-8)")),
       y=expression(paste("Effective Population Size ",N[e]," (x", 10^4 ,")")))
 
#ggsave(smc_bootstrap_plot,filename = "Fig.2A.jpg",width = 8.2,height = 4.3)

climate_data <- read_tsv("data/hpc/demography/nature07158-s2.txt",skip = 14)
cp <- ggplot(climate_data %>% filter(Time>1),aes(x=Time*1e3,y=-Ice_tot)) + 
  geom_line() + 
  scale_x_log10(breaks=c(1e+3,1e+4,1e+5,1e+6), labels=c("1kya","10kya","100kya","1mya"), limits=c(1e+3,1.3e+6),expand = c(0, 0)) + 
  theme_minimal_grid() + labs(x="") +
  ylab("Sea Level")
#ggsave(cp,filename = "Fig.2B.jpg",width = 8.2,height = 2.6)


## this needs to be fixed to plot bootstrap results
bs_ll <- read_tsv("data/hpc/demography/fastsimcoal/bootstrap_param.txt") %>% pivot_longer(-c(MaxEstLhood),names_to = "param",values_to = "value")
  
bs_plot <- bs_ll %>% filter(grepl(param,pattern = "^TDIV")) %>% 
  ggplot(aes(x=param,y=value*5)) +geom_boxplot(aes(fill=param),size=.3,outlier.size = .5) + coord_flip() + theme_cowplot() +
  scale_y_log10(breaks=c(1e+3,1e+4,1e+5,1e+6), labels=c("1kya","10kya","100kya","1mya"), limits=c(1e+3,1.3e+6),expand = c(0, 0)) +
  labs(x="",y="") + scale_fill_brewer(palette="Dark2") + theme(legend.position = "none")

#bs_ll %>% filter(grepl(param, pattern = "^MIG")) %>%  ggplot(aes(x=param,y=value)) +geom_boxplot(aes(fill=param),size=.3,outlier.size = .5)
#ggsave(bs_plot,filename = "Fig.2C.jpg",width = 8.2,height = 3)

library(cowplot)
plot_grid(smc_bootstrap_plot,cp,bs_plot, ncol = 1, align = "hv", axis = "lr", rel_heights = c(0.5,0.25,0.25), labels = c("A","B","C"))
#ggsave("Fig2.jpg",width = 8.2, height = 8.6)

##D
summarise_ld <- function(filename) {
  df <- read_table2(filename) %>% select(-X8) %>% as_tibble %>% 
    mutate(dist = BP_B-BP_A) %>% arrange(dist) %>%
    mutate(distc = cut(dist,seq(min(dist)-1,max(dist)+1,by=200))) %>% 
    group_by(distc) %>% summarise(mean_dist=mean(dist),mean_r2=mean(R2))
  df
}

noffshore_ld <- summarise_ld("data/hpc/popgen/northoffshore.ld") %>% add_column(location="north offshore")
soffshore_ld <- summarise_ld("data/hpc/popgen/southoffshore.ld") %>% add_column(location="south offshore")
inshore_ld <- summarise_ld("data/hpc/popgen/inshore.ld") %>% add_column(location="inshore")


ld_plot <- ggplot(rbind(inshore_ld,noffshore_ld,soffshore_ld),aes(mean_dist/1000,mean_r2)) + 
  geom_point(size = .2,aes(color=location)) + 
  ggsci::scale_color_startrek() +
  theme_test(base_size = 14) + 
  theme(legend.position = "none",
        legend.title = element_blank()) + labs(x="Pairwise distance in Kb", y="LD (" ~ r^2 ~ ")")
ggsave(ld_plot,filename = "Fig.2D.jpg",width = 3,height = 2.5)





## E
read_td <- function(pop) {
  read_tsv(paste0("data/hpc/popgen/",pop,".filtered_windowed.td"), col_names = c("chr","start","end","td")) %>% 
    add_column(pop=pop)
}
td<- rbind(read_td("inshore"),read_td("northoffshore"),read_td("southoffshore"))
td_plot <- ggplot(td,aes(x=pop,y=td,fill=pop)) +
  geom_boxplot(outlier.shape = NA,) + 
  ggsci::scale_fill_startrek() + theme_test(base_size = 14) + 
  scale_y_continuous(limits = c(-3,2)) +
  theme(legend.position = "none") + 
  labs(x="",y="Tajima's D") + 
  scale_x_discrete(labels=c("IN","NO","SO"))
ggsave(td_plot,filename = "Fig.2E.jpg", width = 2.7, height = 2.5)


## F
hbd <- read_tsv("data/hpc/ibdseq/1M_scaffolds.hbd",col_names = c("s1","h1","s2","h2","chr","start","end","LOD")) %>% 
  mutate(loc = locations[str_replace(s1,"_merged","")]) %>% 
  mutate(roh_len = (end-start)/1e6)
source("scripts/color_scheme.R")
hbd_plot <- hbd %>% 
  group_by(s1,loc) %>% 
  summarise(hbd_len = sum(roh_len)) %>% 
  ggplot(aes(x=loc,y=hbd_len)) + 
  geom_boxplot(aes(color=loc),alpha=0.5) +
  geom_point(aes(x=loc,color=loc), position = position_jitter(w=0.1,h=0),size=.6) + 
  xlab("") + ylab("Total length of RoH in Mb") +
  theme_test(base_size = 14) +
  theme(legend.position = "none", text = element_text(size=12)) +
  scale_color_manual(values = myCol) +
  scale_x_discrete(labels=c("IN","NO","SO"))

ggsave(hbd_plot,filename = "Fig.2F.jpg",width = 2.5,height = 2.5)

