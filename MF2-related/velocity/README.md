---

title: RNA velocity analysis for organizer lineages
author: Shaojie Ma
date: Mar 18, 2022

---


#### Prepare the data
- subset to specific subtypes
    - lineages: 
        1. anterior ventral: NKX2-1 & GnRH1+ neurons
        2. rostral patterning center & TAGLN3+ neurons
        3. cortical hem & cajal retzius cells
        4. antihem & IPC TCF7L2
    - script: lin.gen.v2.R & lin.gen.sup.R

- Convert seurat object to h5ad
    - scripts: lin.toh5ad.R & lin.toh5ad.sup.R


#### velocity analysis
- scvelo analysis (stream lines):
    - scripts:
        1. lin.velo.stream.antihem.py
        2. lin.velo.stream.av.py
        3. lin.velo.stream.fgf17.py
        4. lin.velo.stream.hem.py

- scvelo heatmap (gene cascades):
    - pseudotime scripts:
        1. lin.velo.ptime.antihem.py; 
        2. lin.velo.ptime.av.py; 
        3. lin.velo.ptime.fgf17.py; 
        4. lin.velo.ptime.hem.py; 
    - heatmap scripts:
        1. lin.plot.heat.R/sh



