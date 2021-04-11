# Supplementary note

### Sample sequencing and data processing
- [Sampling sites and sequencing data](01.sample_information.md)
- [Variant Calling and quality control](02.quality_control.md)
- [Phasing haplotypes](03.phasing.md)

### Population genetics analyses
- [Population genomics genetic statistics](03.popgen_stats.md)
- [Population structure analysis](04.population_structure.md) and [fineStructure analysis](fs.md)
- [IBD segments](05.ibd_segments.md)
- [Demographic history](05.demographic_history.md)
- [Selection analysis](06.selection_analysis.md)
- [Outgroup analysis](10.outgroup_analysis.md)

### Parallelisation
Some codes provided markdown chunks show how to run command for a single scaffold or a single sample. Most jobs that require to run in HPC were performed using `snakemake` or `parallel` to achieve parallelisation. The actually scripts can be found in [smk](scripts/smk) or [bash](scripts/bash).

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

# Optional download (large files)
wget 'https://cloudstor.aarnet.edu.au/plus/s/Y9Vjkmz3BLL5ogd/download' -O data_large.tgz
tar -zxvf data_large.tgz 
```

Open the project file in RStudio and open the desired file. After installing any required R packages the code should run and produce plots and tables identical to those shown in the web links above.

### Special note
Sometimes, jobs will fail in `snakemake` pipeline, this is common in `shapeit` or `selscan` when there is no SNPs in certain scaffold left after the filtering step and it will raise an error. You can `----keep-going` or `-k` in `snakemake` command line to keep independent jobs going and simply ignore the failed results.
