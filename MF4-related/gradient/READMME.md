---

title: RGC analysis (region gradient)
author: Shaojie Ma
date: August 9th, 2022

---


## Region gradient, Genes enriched in more than one region

#### Type-1. Gene expression are comparable in the enriched regions
- First perform pair-wise region DE comparisons \
deg.region.pairwise.calc.v1.R \
Note here it's pair-wise comparison across regions, while for the region-specific cascades, it's comparison between one region versus all others together. 

- Test region enrichment model (e.g., if a gene is enriched in two regions) \
enrich.model.test.v1.R \

- Extract genes with region-shared patterns & also requiring correlation coefficient \
enrich.mode.extract.reg-share.R \

- Order the genes along pseudotime: \
order.expr.byexp.pairwise.share-reg.trg.R \

- Visualization: \
region.gradients.type-1.cascades.ipynb \


#### Type-2. Gene expression exhibit temporal differences between regions, eg. PENK
region.gradients.type-2.anticor.ipynb \



