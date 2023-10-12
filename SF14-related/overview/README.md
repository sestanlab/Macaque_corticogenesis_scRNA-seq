---

title: Inhibitory neuron (Monkey cortical development)
author: Shaojie Ma
date: July 29, 2022

---


#### Organize InN data from all ages
- script: org.InN.R 


#### Integrate InNs across all ages
- Integrate IPC-InN subtypes across all ages 
    - script: inte.allinn.harmony.v1.R 

- Visualization
    - script: allin.harmony.view


#### Match with adult dataset
- Prepare adult data
    - inte.fetal.adult.pre.adult.R

- Integration 
    - inte.fetal.adult.inn.R/sh 

- Visualization 
    - inte.fetal.adult.view.ipynb


#### Integrate the filtered InNs with Pollen's data. 
- download: The data was downloaded from NCBI GEO, GSE169122. 

- organize data: inte.pollen.org.data.v2.R 

- integration: 
    - notes: Both Seurat and harmony batch correction methods were tried. Seurat is slow but has better performance
    - scripts: inte.pollen.seurat.R/sh & inte.pollen.harmony.v3.R/sh 

- Visualization: inte.pollen.view.ipynb













