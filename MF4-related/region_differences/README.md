---

title: RGC analysis (region differences)
author: Shaojie Ma
date: August 9th, 2022

---



### Region-specific cascade of gene expression

#### Find region DEGs
- calculation: 
    - script: cascade.deg.region.calc.v1.R 

- extract DEGs: 
    - script: cascade.deg.region.extract.v1.R 

- Visualize the number of DEGs along the pseudotime
    - prepare the data for matlab visualization:
        - prepare the seurat object: sliding.degs.calc.org.prep.R, sliding.degs.calc.trg.prep.R
        - get the nDEGs table for matlab visualization: sliding.degs.calc.org.R, sliding.degs.calc.trg.R

    - matlab visualization (3D ribbon plot): 
        - matlab_div_nDEGs_oRG_v2.m
        - matlab_div_nDEGs_tRG_v2.m
        - source functions: ./matlab_custom_scripts/ribboncoloredZ.m


#### Calculate average expression & smooth
- Average expression by bins: 
    - script: cascade.avg.region.calc.v1.R

- Smooth gene expression along pseudotime: 
    - script: cascade.avg.smooth.region.calc.v1.R


#### Order genes along the trajectory
- Order genes based on their expression (smoothed expression) in a given region
cascade.order.expr.byexp.v4.sh \
cascade.order.expr.byexp.org.v4.R \
cascade.order.expr.byexp.trg.v4.R \


#### Visualization
- NESC-tRG MEST lineage: cascades.view.tRG.ipynb
    - show all genes
    - show genes in selected signaling pathways

- oRG lineage: cascades.view.oRG.ipynb
    - show all genes
    - show genes in selected signaling pathways


#### oRG makers utilized by early frontal RGCs
- Evaluate whether oRG markers can distinguish early FC versus OC RGC subtypes 
    - script: enrich.org.marker.frontal.early.nsc.v2.ipynb


#### DEGs shared between early and late RGCs
- scripts: degs.shared.early.late.rgcs.ipynb







