---
title: Bulk tissue analysis
author: Shaojie Ma
date: Nov 23, 2022
---


#### Organize the data needed for the analysis
- Bulk tissue
    - Data are shared from Mingfeng
    - script: organize.data.bulk.ipynb

- single cell
    - subset from the full dataset
    - script: organize.data.vivo.v2.ipynb


#### DeSeq2 of the bulk tissue RNA-seq data
- Concept
    - in vitro data were divided to three stages (early, middle, late) and we wanted to check how many DEGs are present in each stage and how they are overlapped 
    - script: de.bulk.v2.ipynb. In this scirpt, the overlaps of region-specific DEGs across the three stages were visualized.


#### Connect in vitro and in vivo results
- Concept: We categorized two types of genes using the in vitro system
    - type1: validate if the in vivo region-specific features of NSCs are maintained in vitro
    - type2: genes that are only upregulated in the late stage (more mature ExNs)
    - analysis: then we check for each type genes, whether we can see similar patterns in vivo

- script: bulk.degs.nsc.vitro.vivo.ipynb
    - in vivo: check if any of the in vitro region-specific genes are also displaying the same region-specific expression in vivo NSCs
    
- script: bulk.degs.type2.invivo.v3.ipynb
    - in vitro: check if any of the type2 genes are upregulated in the same region in differentiated neurons. 


