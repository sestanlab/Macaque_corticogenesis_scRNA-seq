---
title: TF-target network for the brain organizers
author: Shaojie Ma
date: Feb 15, 2022
---


## Source functions
- pwm_codes.R: PWM processing functions
- motif.fun.R: motif enrichment functions
- net.fun.v2.R: network visualization functions


## Motif enrichment analysis
- Build motif enrichment background for PWMEnrich
	- motif.build.bg.R: 
		1. use the expressed genes to extract all the related PWMs
		2. extract the promoter regions of the genome
		3. build the background


## Prepare for the network analysis
- calculate average expression along UMAP axis
	- net.actual.avg.R: the actual average expression
	- net.perm.avg.R: average expression of the permutated data. The permutation helps to decide how robust the co-expression are.


## Network calculation and visualization
net.signaling.v2.ipynb


