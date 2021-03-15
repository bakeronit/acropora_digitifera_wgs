## Supplementary information

### Sample sequencing and data processing
- [Sampling sites and sequencing data](01.sample_information.md)
- [Variant Calling and quality control](02.quality_control.md)
- [Phasing haplotypes]

### Population genetics analyses

- [Population structure analysis](04.population_structure.md)


All of the sections above are provided as processed markdown files. Clicking the link should display a web readable page with text, a few select commands and plots and tables. The code used to generate these pages is provided in the corresponding .Rmd file. If you would like to run the code in these files yourself you will need to;

Checkout this repository

```bash
git clone https://github.com/bakeronit/acropora_digitifera_wgs.git
```

Download additional data not hosted on github due to size constraints

```bash
cd acropora_digitifera_wgs
wget https://cloudstor.aarnet.edu.au/plus/s/j9hwbpnYKbNUWTL/download -O data.tgz
tar -zxvf data.tgz 
```

Open the project file in RStudio and open the desired file. After installing any required R packages the code should run and produce plots and tables identical to those shown in the web links above.
