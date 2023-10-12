---

title: Regional differences of ExN lineage
author: Shaojie Ma
date: Aug 31, 2022

---


#### Number of DEGs along pseudotime
- prepare data: ndegs.prep.data.R

- calculation: ndegs.calc.deep.R & ndegs.calc.upper.R

- visualization in matlab: matlab_IE_div_IT.m & matlab_IE_div_nonIT.m


#### Augur measuring regional differences of the shared cell subtypes
- prepare the data
    - organize data: augur.pre.data.R
    - prepare the genes (HVGs): augur.pre.hvg.R

- run Augur
    - deep layer: augur.run.deep.R
    - upper layer: augur.run.upper.R

- visualization
    - augur.view.ipynb


#### Hierachical clustering of the same subtype across regions
- prepare data
    - organize data: hc-reg.predata.R
    - HVGs: hc-reg.prehvg.R

- permutated data
    - permutated AVGs: hc-reg.perm.gen.R/sh

- Visualization
    - script: hc_hierarchical-clustering.region.visualization.ipynb


#### Gene expression cascades across regions
- prepare data
    - pseudobulk avgs along pseudotime: prep.avgs.up.v1.R & prep.avgs.deep.v1.R 
    - degs: prep.degs.deep.v1.R & prep.degs.up.v1.R
- visualization:
    - upper: region-specific.gene.cascades.upper.ipynb
    - deep: region-specific.gene.cascades.deep.ipynb


#### Expression gradients across regions
- upper: expression.gradient.upper.R

- deep: expression.gradient.deep.R


#### Correlation of region-specificity among ExN subtypes
- prepare DEGs: 
    - scripts: mars.calc.upper.R; mars.calc.deep.R

- calculate correlation: 
    - scripts: correlate.degs.deep-upper.ipynb


#### DEGs shared between RGC subtypes and ExN subtypes
- script: DEGs.share.RGC.ExN.ipynb




