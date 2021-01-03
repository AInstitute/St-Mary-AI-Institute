---
title: Machine Learning For Tracing the Transmission of SARS-Cov-2
date: 2020-12-16T09:51:12.000+00:00
image_webp: "/images/105-72.jpg"
image: "/images/haplotype_ncov2019_from_ngdc.svg"
author: Jack Zhang
description: Data Visualization

---
# Spatial genetic structure of SARS-CoV-2

## An interactive map of global spatial genetic structure of SARS-CoV-2

To investigate the spatial genetic variations across these 62 regions, we conducted principal component analysis on the genotypes of these virus samples. We first did the PCA on the bi-allelic markers (SNPs). We exacted 2, 551 bi-allelic markers (SNPs) and conducted the PCA with minor allele frequency (MAF) >0.05 and LD at r2 =0.8 using the SNPrelate package (Zheng, et al. 2012). We also carried out the PCA on all SNPs in order to capture a thorough understanding of genome-wide single-nucleotide polymorphism variations. As this virus was reported to have a high mutation rate (Phan 2020; Shen, et al. 2020), to construct a more comprehensive genetic map across space, we further carried out the principal component analysis on all the variant sites (CNRs, the copy number of references), including bi-allelic SNPs, multi-allelic SNPs, indels and structural variants. To see the genetic structure and the genetic relation of these viruses from some heavily impacted countries as well as the early reported countries, we also presented the genetic variation map for China, America, Australia, and Canada. We provide interactive spatial genetic maps for two types of markers (SNPs and CNRs) for users to look into the local and global structure of this virus.

### Global spatial genetic structure of SARS-CoV-2 represented by SNPs

![](/images/fig1-1-1.png)

Fig. 1. Spatial genetic structure of SARS-CoV-2 based on whole genome wide SNPs

Readers can use their mouses to zoom in and zoom out to check the samples. Double click the legend will hidde the corresponding points from the plots.

### Global spatial genetic structure of SARS-CoV-2 represented by copy number of variations

![](/images/fig2-1-1.png)

Fig. 2. Spatial genetic structure of SARS-CoV-2 based on whole genome wide copy number of variations.

#### Reference

###### Data source: CNCB- 2019 Novel Coronavirus Resource (2019nCoVR).

###### Source:[https://xinghuq.github.io/SARS-COV-2/Spatial_structure.html](https://xinghuq.github.io/SARS-COV-2/Spatial_structure.html "https://xinghuq.github.io/SARS-COV-2/Spatial_structure.html")