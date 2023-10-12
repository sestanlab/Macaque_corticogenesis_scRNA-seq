---

title: RGC analysis (overview)
author: Shaojie Ma
date: August 9th, 2022

---


#### Markers
- some selected script: 
	- script: marker.view.ipynb
	    - visualize known major markers
	    - visualize the markers for the frontal-specific subtype
    - metadata for the clusters: plot.frontal-spec.cluster.mar.meta.ipynb

- shared gene cascades across regions
    - prepare data:
        - average expression: pseudobulk.avg.pseudotime.R
        - calculate markers: shared.progression.marker.calc.R
        - order genes based on impulsefit: shared.progression.marker.urd.order.R
    - visualization: 
        - script: shared.progression.marker.cascades.ipynb


#### Data comparisons with public data
- script: match.public.rgc.ipynb
    - Comparison with Bhaduri et al., 2021 Nature
    - Comparison with nowakowski et al., 2018 Science





