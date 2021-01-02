---
title: Using machine learning to trace the recent transmission of SARS-Cov-2
date: 2020-11-03T09:51:12+00:00
author: Mark Jodan
image_webp: "/images/haplotype_ncov2019_from_ngdc.svg"
image: "/images/da_example.png"
description: Data visualization

---
## Background

Successive transmission of SARS-CoV-2 has caused a devastating influence on human health and the global economy. Previous studies have reported the origin, genetic mutation, phylogenetic network, etc, of SARS-CoV-2, since the first report of its outbreak in Wuhan. However, spatial genomic analysis, from a population genomic perspective, has not entirely demonstrated its power in tracing the origin and transmission, as well as the evolution of SARS-CoV-2 so far.

## Methods

We identified genome-wide variant sites from 3,736 complete genomes of SARS-CoV-2 isolates from 62 geographic regions. To trace the origin and migration of SARS-CoV-2 on a local and regional basis, we carried out spatial genetic analysis of whole-genome variants, including single-nucleotide polymorphisms (SNPs) and copy number of references (CNRs). Further, we conducted the spatial variance-based genome scan on these genome-wide variants using pcadapt to detect the genome regions under selection or involved in host adaptation.

 Fig. 1. SARS-Cov-2 sampling regions. The abbreviations correspond to the ISO country (region) code.

## Findings

Spatial genetic variations indicate that the viruses from 62 regions presented a significant genetic structure with several clusters (four clusters in SNPs, and three clusters in CNRs) showing strong genetic differentiation. Notably, whole-genome SNP analysis revealed a unique cluster in America, Canada, and Australia isolates, while this cluster was not found in other geographic regions, indicating the possibility of isolation due to distance or through selection. Analysis of genome-wide copy number of references (CNRs) uncovered the nuanced spatial genetic variations of SARS-CoV-2 associated with migration and transmission. Genome scan results showed that SARS-CoV-2 has an unprecedented rapid evolution with numbers of loci that were potentially under selection. Remarkably, S, E, and M genes showed evidence of strong selection than other genes among genome-wide variant sites.

## Interpretation

Linear spatial cline between neighbour countries indicates the high transmission possibility caused by frequent migration. Therefore, it is effective to constrain human migration (quarantine) to prevent recurrent spread. It is still too early to conclude the origin of the SARS-CoV-2. However, it is clear that currently available data suggest that viruses in China acted as genetic chains connecting to the divergent sub-isolates, even to the origin. The accumulation of more comprehensive data from various isolates, including these from different animals, will dig out the origin in the near future through our inspection approach. Our results also suggested that host selection on genetic variations increases the challenge of medical development and treatment. But again, it is important to consider these variants to take specific or personalized clinical strategies to treat patients infected by different sub-strains or through poly-targeted drugs.

## Datasets

The genome sequences and genome variations were obtained from CNCB - 2019 Novel Coronavirus Resource (2019nCoVR) (Zhao, et al. 2020)( [https://bigd.big.ac.cn/ncov?lang=en](https://bigd.big.ac.cn/ncov?lang=en "https://bigd.big.ac.cn/ncov?lang=en")), where complete genome sequences used for analyses within this resource were obtained from the CNGBdb (Wang, et al. 2019), GenBank (Benson, et al. 2011), GISAID (Shu and McCauley 2017), GWH (Center 2020) and NMDC databases ([https://microbiomedata.org/](https://bigd.big.ac.cn/ncov?lang=en "https://bigd.big.ac.cn/ncov?lang=en")). This data includes 3,736 complete genomes of SARS-CoV-2 strains isolated from human individuals from 62 regions (Fig. 1), including countries, administrative regions, and one cruise ship. We identified 3, 218 variant sites in total, including bi-allelic SNPs, multi-allelic SNPs, indels, and structural variants.

CITATION: Qin. X, 2020, Spatial genetic structure of SARS-CoV-2.SARS-COV-2.[https://github.com/xinghuq/SARS-COV-2.](https://bigd.big.ac.cn/ncov?lang=en "https://bigd.big.ac.cn/ncov?lang=en")