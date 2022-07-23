## Population genetics of *Acropora digitifera* in North-Western Australia using whole genome sequencing

This repository contains complete computational methods including source code, processed data outputs and explanatory notes.  It is an accompaniment to the paper;

> Publication link TBA


### Table of Contents

**Sample sequencing and data processing**
- [Sampling sites and sequencing data](01.sample_information.md)
- [Variant Calling and quality control](02.quality_control.md)
- [Phasing haplotypes](03.phasing.md)
- [Building pseudo-chromosomes](x30.ragtag_scaffolding.md)

**Population genetics of the coral host**
- [Basic population genetic statistics](04.popgen_stats.md)
- [Population structure analysis](05.population_structure.md)
- [FineStructure analysis](x20.finestructure.md)
- [Outgroup analysis using Japanese A. digitifera](x10.outgroup_analysis.md)
- [IBD and HBD segments](06.ibd_hbd.md)
- [Demographic history based on SMC](07.demographic_history.md)
- [Demographic model fitting using fastsimcoal](22a.fastsimcoal_fitting.md)
- [Demographic model validation by simulation](22b.fastsimcoal_sim.md)
- [Demographic model without a bottleneck](24.nobottle_sim.md)
- [All .est and .tpl files used for alternative demographic models](data/models/)

**Symbiodiniaceae profiles**
- [Dissecting the symbiont reads](23.symbionts.md)

**Signatures of selection in the coral host**
- Selection analysis based on [EHH](08.ehh_stats.md) and [allele frequency](12.pbs.md)
- [Comparison between EHH and allele-frequency-based approaches](14.ehh_pbs.md)
- [Identification of the top 1% of regions putatively under selection](10.identify_selective_genomic_windows.md)
- [GO Enrichment of genes under selection](15.GO_enrichment.md)
- [Estimating the TMRCA for selective sweeps](17.dating_the_selection.md)
- [Phylogenetic tree of haem peroxidases](27.peroxidases.md)


**Miscellaneous**
- [Functional annotation of genes](09.annotate_genes.md)
- [Checking the mislocated sample using radseq data](18.radseq_check.md)
- [Estimating the age of individual alleles](25.geva.md)


### How to use this repository
All of the sections above are provided as processed markdown files. Clicking the link should display a web readable page with text, a few select commands and plots and tables. The code used to generate these pages is provided in the corresponding .Rmd file. If you would like to run the code in these files yourself you will need to;

**Checkout this repository**

```bash
git clone https://github.com/bakeronit/acropora_digitifera_wgs.git
```

Download additional data not hosted on github due to size constraints

```bash
cd acropora_digitifera_wgs
wget 'https://cloudstor.aarnet.edu.au/plus/s/6GcW7OUebjrgASm/download' -O data_essential.tgz
tar -zxvf data_essential.tgz 
```


*Optional:* download (large files).  Includes the full phased and unphased SNP sets

```bash
wget 'https://cloudstor.aarnet.edu.au/plus/s/Y9Vjkmz3BLL5ogd/download' -O data_large.tgz
tar -zxvf data_large.tgz 
```

Open the project file in RStudio and open the desired file. After installing any required R packages the code should run and produce plots and tables identical to those shown in the web links above.

### Code Portability
In general the code provided here is for illustration purposes only.  In many cases the commands we provide should work as-is on your system.  In some cases however, we have provided scripts that include strategies for parallelisation (eg via gnu parallel, or using a PBS job management system).  These are highly platform specific and you will need to modify them to run on your own system.

