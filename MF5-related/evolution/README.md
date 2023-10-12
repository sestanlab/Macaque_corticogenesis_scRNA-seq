---
title: Evolutionary comparisons between human and macaque data
author: Shaojie Ma
date: Nov 1, 2022
---


#### Organize the data
- Create the seurat objects for human and macaque
    - script: organize.data.R
- downampling
    - script: downsample.data.R
    - concept: both human and macaque data were downsampled to the same number of UMIs (1500) to avoid sequencing bias
    - source functions: ds.fun.R


#### Overall transcriptomic similarity between human and macaque subtypes
- script: tp.simi.human.macaque.ipynb
- concept: use shared HVGs & pearson correlation to denote the overall transcriptomic similarity


#### Species comparisons (shared and species-specific region-specific features)
- script: compare.human.macaque.ipynb
- source functions: vis.fun.R
- concept: use the downsampled data to do the comparisons. Find the top genes shared or enriched in a given species



