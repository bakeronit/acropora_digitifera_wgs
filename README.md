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
- [IBD and HBD segments](06.ibd_hbd.md)
- [Demographic history based on SMC](07.demographic_history.md)
- [Demographic model fitting using fastsimcoal](22a.fastsimcoal_fitting.md)
- [Demographic model validation by simulation](22b.fastsimcoal_sim.md)
- Selection analysis based on [EHH](08.ehh_stats.md) and [allele frequency](12.pcangsd_selection.md)
- [GO Enrichment of genes under selection](11.GO_enrichment.md)
- [GO Enrichment of genes under selection based on Interproscan GO term assignments](x40.GO_ipr_enrichment.md)
- [Estimating the TMRCA for selective sweeps](17.dating_the_selection.md)

**The symbiodiniacea profies**
- [Dissecting the symbiont reads](23.symbionts.md)

**Things waiting to be ordered**
- [Gene function annotation](09.annotate_genes.md)
- [EHH windowed scan, how I get the top 1% genomics region from iHS, XP-xxx](10.identify_selective_genomic_windows.md)
- [Pcangsd_popgen_stats, Fst/TajimasD using ANGSD](13.popgen_stats_angsd.md)
- [Trend that EHH candidate region with high pbs and high tajimasD](14.ehh_pbs_pcangsd.md)
- [Checkind the mislocated sample using radseq data](18.radseq_check.md)
- [Combined EHH scan candidate in inshore](20.inshore_candidate_genes.md)
- [Jia's old GO enrichment analysis](21.functional_enrichment.md)
- [outgroup analysis, check using JP a.digitifera](x10.outgroup_analysis.md)
- [fineStructure analysis](x20.finestructure.md)
- [Ira's smc plot with mutation rate range](x50.smc_files.md)
- [host mitogenome haplotype](24.host_mitogenome.md)

**Unwanted Rmd or not rendered Rmd**
- [19.ldheatmap.Rmd](19.ldheatmap.Rmd): delete it?
- [16.symbiont_profiles.Rmd](16.symbiont_profiles.Rmd): replaced by 23.symbiont.Rmd
- [15.pcangsd_neutral.Rmd](15.pcangsd_neutral.Rmd)

### How to use this repository
All of the sections above are provided as processed markdown files. Clicking the link should display a web readable page with text, a few select commands and plots and tables. The code used to generate these pages is provided in the corresponding .Rmd file. If you would like to run the code in these files yourself you will need to;

**Checkout this repository**

```bash
git clone https://github.com/bakeronit/acropora_digitifera_wgs.git
```

Download additional data not hosted on github due to size constraints

```bash
cd acropora_digitifera_wgs
wget 'https://cloudstor.aarnet.edu.au/plus/s/BKYCAsujh2sjLdz/download' -O data.tgz
tar -zxvf data.tgz 

*Optional download (large files)*
wget 'https://cloudstor.aarnet.edu.au/plus/s/Y9Vjkmz3BLL5ogd/download' -O data_large.tgz
tar -zxvf data_large.tgz 
```

Open the project file in RStudio and open the desired file. After installing any required R packages the code should run and produce plots and tables identical to those shown in the web links above.

### Parallelisation in HPC
Some codes provided in the markdown bash chunks show how to run command for a single scaffold or a single sample. Most jobs that require to run in HPC were performed using `snakemake` or GNU `parallel` to achieve parallelisation. The actually scripts can be found in [smk](scripts/smk) or [bash](scripts/bash). Sometimes, jobs will fail in `snakemake` pipeline, this is common in `shapeit` or `selscan` when there is no SNPs in certain scaffold left after the filtering step and it will raise an error. You can `--keep-going` or `-k` in `snakemake` command line to keep independent jobs going and simply ignore the failed results.
