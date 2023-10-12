---

title: Inhibitory neuron (Monkey cortical development)
author: Shaojie Ma
date: July 29, 2022

---



#### Augur analysis to discriminate regional differences
- Augur prepare data 
run.augur.pre.data.R 

- Run Augur for each region-pair in each age \
run.augur.R/sh \

- Plot Augur results
run.augur.plot.R \



#### Do differential gene expression
- Calculate avg & DEGs: deg.avg.calc.v1.R 

- script: deg.region.correlation.ipynb
    - Show region expression gradient
    - Correlation of region enrichment (for both MGE/CGE InN)




