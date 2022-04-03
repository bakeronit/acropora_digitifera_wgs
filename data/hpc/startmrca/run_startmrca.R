library(tidyverse)
library(rehh)
library(startmrca)


args <- commandArgs(trailingOnly = TRUE)

position <- args[1]
derived_allele <- as.integer(args[2]) #0 or 1 if iHS is positive it will be 1, otherwise 0
vcf_file <- args[3]
sweep_name <- args[4]
pop <- args[5]
rep <- as.integer(args[6])

cat("Reading ",vcf_file," \n")


hh <- data2haplohh(hap_file = vcf_file,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

marker_index <- which(positions(hh)==position)

# Create lists of individuals for ancestral and derived states
#


derived_haps <- which(haplo(hh)[,marker_index]==derived_allele) %>% names()
anc_haps <- which(haplo(hh)[,marker_index]!=derived_allele) %>% names()

cat(marker_index," ",position,"\n")

# derived_indv <- derived_haps %>% str_remove("_[0-9]$")
# anc_indv <- anc_haps %>% str_remove("_[0-9]$")

#derived_only_indv <- setdiff(derived_indv,anc_indv) %>% as.data.frame()
#anc_only_indv <- setdiff(anc_indv,derived_indv) %>% as.data.frame()


derived_path <- paste(pop,"_indvs.txt",sep="")
anc_path <- paste(pop,"_ref_indvs.txt",sep="")



# Run startmrca

# For nsel and nanc we want at least 10 in each and at most 40
if ( (length(derived_haps) < 6) | (length(anc_haps) < 6) ){
    stop("Not enough haplotypes in ancestral or derived to perform sampling")
}

nsel = length(derived_haps)
nanc = length(anc_haps)

cat("Number of selected haps: ",nsel," and ancestral: ",nanc,"\n")


#for (rep in 1:3) {
    tmrca_results <- run.startmrca(vcf.file=vcf_file, 
                        rec.file=NULL, 
                        mut.rate = 1.2e-8,rec.rate= 3.2e-8, 
                        nsel = nsel, nanc = nanc, 
                        chain.length = 10000, 
                        proposal.sd = 20, 
                        nanc.post = 100, 
                        sample.ids = derived_path, 
                        refsample.ids = anc_path,
                        pos = position, 
                        sel.allele = derived_allele, bed.file = NULL, upper.t.limit = 10000, 
                        output.file = paste("results/",sweep_name,"_mcmc_list.RDATA",sep=""))

    write_rds(tmrca_results,paste("results/",sweep_name,"_",rep,".rds",sep=""))
#}


