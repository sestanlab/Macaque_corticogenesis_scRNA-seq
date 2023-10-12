---

title: Disease analysis (RGC progression)
author: Shaojie Ma
date: July 29, 2022

---



#### RGC Pseudobulk AUCell enrichment
- prepare data
    - data: org.rgc.data.R
    - disease genes: org.disease.genes.v2.R

- Generate pseudobulk average expression & AUC enrichment
calc.pseudobulk.avg.rgc.v2.R/sh

- Subset the data & generate pseudobulk AVG & AUC enrichment
calc.pseudobulk.avg.subset.rgc.v2.R

- Visualization: rgc.pseudobulk.enrichment.view
    - data: visualize both the raw data and subset data
    - plot: enrichment & meta: age proportions of the RGC subtypes


#### Intersect region-specific genes with disease genes
- prepare data
    - scripts: reg-deg.disease.gene.R & reg-deg.disease.mar.R

- source function: reg-deg.fun.R & pie.fun.R

- visualization: region.deg.markers.view.ipynb


