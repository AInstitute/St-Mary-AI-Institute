---
title: " Ecological and Evolutionary Inference using Supervised Learning"
date: 2020-08-06T06:52:36+00:00
image_webp: "/images/logo.png"
image: "/images/backgrounds/st-andrews2.jpg"
author: Jack Zhang
description: DA2

---
# Community structure inference using discriminant analysis

#### Xinghu Qin –School of Biology, University of St Andrews

#### 03-12-2020

Source: [https://xinghuq.github.io/DA/articles/Microbiome.html](https://xinghuq.github.io/DA/articles/Microbiome.html "https://xinghuq.github.io/DA/articles/Microbiome.html")

## Community structure inference using DA

We will be using data from Wagner et al. (2016), which studies the effects of plant age, genotype, and environment on the bacterial microbiome of a perennial herb, Boechera stricta, in the mustard family. The raw data of Wagner et al. (2016) is available on [dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.g60r3).

I have stored the data into the package and users can load the data from package.

    library(DA)
    
    f <- system.file('extdata',package='DA')
    infile <- file.path(f, "microbiome.Rdata")
    load(infile)
    # Warning: namespace 'taxa' is not available and has been replaced
    # by .GlobalEnv when processing object 'obj'
    micro_abund=as.data.frame(t(obj$data$otu_rarefied[, sample_data$SampleID]))
    
    sample_data$Site=factor(sample_data$Site,levels = unique(sample_data$Site))

## Principal Component Analysis (PCA)

We first use PCA to analyse the bacterial community from the plants at different sites.

    ### PCA
    Microbiome_pc=prcomp(micro_abund,scale. = TRUE)
    
    #plot the data projection on the components
    library(plotly)
    # Loading required package: ggplot2
    # 
    # Attaching package: 'plotly'
    # The following object is masked from 'package:ggplot2':
    # 
    #     last_plot
    # The following object is masked from 'package:stats':
    # 
    #     filter
    # The following object is masked from 'package:graphics':
    # 
    #     layout
       cols=rainbow(length(unique(sample_data$Site)))
       p0 <- plot_ly(as.data.frame(Microbiome_pc$x), x =Microbiome_pc$x[,1], y =Microbiome_pc$x[,2], color = sample_data$Site,colors=cols[sample_data$Site],symbol = sample_data$Site,symbols = 1:3L) %>% 
         add_markers() %>%
         layout(scene = list(xaxis = list(title = 'PC1'),
                             yaxis = list(title = 'PC2')))
    
      
    p0

−20−100102030−10010203040

Fig.1 PCA plot of bacterial microbiome community at different sites.

## Discriminant Analysis of Principal Components (DAPC)

DAPC has widely used in ecology and evolution. Using DAPC to display the community structure of bacterial microbiomes is rarely exploited.

    library(adegenet)
    # Loading required package: ade4
    # Registered S3 method overwritten by 'spdep':
    #   method   from
    #   plot.mst ape
    # Registered S3 methods overwritten by 'vegan':
    #   method      from
    #   plot.rda    klaR
    #   predict.rda klaR
    #   print.rda   klaR
    # 
    #    /// adegenet 2.1.1 is loaded ////////////
    # 
    #    > overview: '?adegenet'
    #    > tutorials/doc/questions: 'adegenetWeb()' 
    #    > bug reports/feature requests: adegenetIssues()
    sample_data$Site=factor(sample_data$Site,levels = unique(sample_data$Site))
    ###DAPC
    Microbiome_dapc=dapc(micro_abund,grp=sample_data$Site,n.pca=10, n.da=3)
    
    #plot the data projection on the components
    library(plotly)
       cols=rainbow(length(unique(sample_data$Site)))
       p1 <- plot_ly(as.data.frame(Microbiome_dapc$ind.coord), x =Microbiome_dapc$ind.coord[,1], y =Microbiome_dapc$ind.coord[,2], color = sample_data$Site,colors=cols[sample_data$Site],symbol = sample_data$Site,symbols = 1:3L) %>% 
         add_markers() %>%
         layout(scene = list(xaxis = list(title = 'DAPC1'),
                             yaxis = list(title = 'DAPC2')))
    
       p1

−20246−4−202468

Fig.2 DAPC plot of bacterial microbiome community at different sites.

This is an interactive plot that allows you to point the data values and display the value as you wish.

## Discriminant Analysis of Kernel Principal Components (DAKPC)

Compared to DAPC, discriminant analysis of kernel principal components (DAKPC) uses the non-liner kernal technique. The kernel principal component analysis is emolyed to incorporate the non-linear relationship between sites and samples in DAKPC. Below is the implementation of DAKPC.

    Microbiome_ldakpc=LDAKPC(micro_abund,sample_data$Site,n.pc=10)
    # Loading required package: kernlab
    # 
    # Attaching package: 'kernlab'
    # The following object is masked from 'package:ggplot2':
    # 
    #     alpha
    
     cols=rainbow(length(unique(sample_data$Site)))
       p2 <- plot_ly(as.data.frame(Microbiome_ldakpc$LDs), x =Microbiome_ldakpc$LDs[,1], y =Microbiome_ldakpc$LDs[,2], color = sample_data$Site,colors=cols[sample_data$Site],symbol = sample_data$Site,symbols = 1:3L) %>% 
         add_markers() %>%
         layout(scene = list(xaxis = list(title = 'LDAKPC1'),
                             yaxis = list(title = 'LDAKPC2')))
    p2

−4−2024−6−4−202468

Fig.3 LDAKPC plot of bacterial microbiome community at different sites.

LDAKPC has the similar result with DAPC.

## Local Fisher Discriminant Analysis (LFDA)

As we mentioned in previous example, LFDA can discriminate the multimodal data while LDA can not. LFDA is an upgraded version of LDA preserving within group variance.

    Microbiome_lfda=LFDA(micro_abund,sample_data$Site,r=3,tol=1E-3)
    # Loading required package: lfda
    # Loading required package: klaR
    # Loading required package: MASS
    # 
    # Attaching package: 'MASS'
    # The following object is masked from 'package:plotly':
    # 
    #     select
    
    cols=rainbow(length(unique(sample_data$Site)))
    p3 <- plot_ly(as.data.frame(Microbiome_lfda$Z), x =Microbiome_lfda$Z[,1], y =Microbiome_lfda$Z[,2], color = sample_data$Site,colors=cols[sample_data$Site],symbol = sample_data$Site,symbols = 1:3L) %>% 
         add_markers() %>%
         layout(scene = list(xaxis = list(title = 'LFDA1'),
                             yaxis = list(title = 'LFDA2')))
    p3

−9−8−7−6−5−4−3−2−1−3−2−1012

Fig.4 LFDA plot of bacterial microbiome community at different sites.

## Local Fisher Discriminant Analysis of Kernel Principal Components (LFDAKPC)

Replacing LFDA for discriminant analysis on the basis of LDAKPC we will get LFDAKPC, Local (Fisher) Discriminant Analysis of Kernel Principal Components (LFDAKPC). Below is the implementation of LFDAKPC.

    Microbiome_lfdakpc=LFDAKPC(micro_abund,sample_data$Site,kernel.name="polydot",kpar = list(degree = 1, scale = 1, offset = 1),n.pc=10,tol=1E-30)
    
    cols=rainbow(length(unique(sample_data$Site)))
    p4 <- plot_ly(as.data.frame(Microbiome_lfdakpc$LDs), x =Microbiome_lfdakpc$LDs[,1], y =Microbiome_lfdakpc$LDs[,2], color = sample_data$Site,colors=cols[sample_data$Site],symbol = sample_data$Site,symbols = 1:3L) %>% 
         add_markers() %>%
         layout(scene = list(xaxis = list(title = 'LFDAKPC1'),
                             yaxis = list(title = 'LFDAKPC2')))
    p4

−5k05k10k15k20k−5000−4000−3000−2000−100001000200030004000

Fig.5 LFDAKPC plot of bacterial microbiome community at different sites.

The LFDAKPC also produces the similar results as LDAKPC and DAPC.

## Kernel Local Fisher Discriminant Analysis (KLFDA)

Kernel local (Fisher) discriminant analysis (KLFDA) capatures the non-linear relationships between samples and also considers the within group multimodality.

    ### Note the kernel matrix used here is from kernelab library, while the kernel techniques in kernelab is different, it uses inverse kernel distance.
    
    Microbiome_klfda=klfda_1(as.matrix(micro_abund),as.matrix(sample_data$Site),kernel=kernlab::rbfdot(sigma = 0.00001),r=3,tol=1E-90,prior = NULL)
    # Loading required package: WMDB
    
    cols=rainbow(length(unique(sample_data$Site)))
    p5 <- plot_ly(as.data.frame(Microbiome_klfda$Z), x =Microbiome_klfda$Z[,1], y =Microbiome_klfda$Z[,2], color = sample_data$Site,colors=cols[sample_data$Site],symbol = sample_data$Site,symbols = 1:3L) %>% 
         add_markers() %>%
         layout(scene = list(xaxis = list(title = 'LDA1'),
                             yaxis = list(title = 'LDA2')))
    p5

−10123−0.6−0.4−0.200.20.40.60.81

Fig.6 KLFDA plot of bacterial microbiome community at different sites.

KLFDA clearly presents the bacterial aggregates between different sites.

All the above methods show the same global structure for bacterial microbiome community at different sites. But LFDA and KLFDA do a better job in discriminating the communities.

## Species assignment using Kernel Local Fisher Discriminant Analysis (KLFDA)

Kernel local (Fisher) discriminant analysis (KLFDA) does the best for analysis of community structurte. Now, we plot the Microbiome individual membership representing the posterior possibilities of species as the community structure.

    
    library(adegenet)
    ## asignment plot
    compoplot(as.matrix(Microbiome_klfda$bayes_assigment$posterior),show.lab = TRUE, posi=list(x=5,y=-0.01),txt.leg = unique(sample_data$Site))

![&nbsp;](https://xinghuq.github.io/DA/articles/Microbiome_files/figure-html/fig7-1.png =816x)

Fig. 7 The community structure of bacterial microbiome community at different sites (individual assignment)

# References

Laloë, D., Jombart, T., Dufour, A.-B. & Moazami-Goudarzi, K. (2007). Consensus genetic structuring and typological value of markers using multiple co-inertia analysis. Genetics Selection Evolution, 39, 545.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24, 1403-1405. Sugiyama, M (2007).Dimensionality reduction of multimodal labeled data by local Fisher discriminant analysis. Journal of Machine Learning Research, vol.8, 1027-1061.

Sugiyama, M (2006). Local Fisher discriminant analysis for supervised dimensionality reduction. In W. W. Cohen and A. Moore (Eds.), Proceedings of 23rd International Conference on Machine Learning (ICML2006), 905-912.

Original Matlab Implementation: [http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA](http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA "http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA")

Tang, Y., & Li, W. (2019). lfda: Local Fisher Discriminant Analysis inR. Journal of Open Source Software, 4(39), 1572.

Moore, A. W. (2004). Naive Bayes Classifiers. In School of Computer Science. Carnegie Mellon University.

Pierre Enel (2020). Kernel Fisher Discriminant Analysis ([https://www.github.com/p-enel/MatlabKFDA](http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA "http://www.ms.k.u-tokyo.ac.jp/software.html#LFDA")), GitHub. Retrieved March 30, 2020.

Karatzoglou, A., Smola, A., Hornik, K., & Zeileis, A. (2004). kernlab-an S4 package for kernel methods in R. Journal of statistical software, 11(9), 1-20.

Bingpei Wu, 2012, WMDB 1.0: Discriminant Analysis Methods by Weight Mahalanobis Distance and bayes.

Ito, Y., Srinivasan, C., & Izumi, H. (2006, September). Discriminant analysis by a neural network with Mahalanobis distance. In International Conference on Artificial Neural Networks (pp. 350-360). Springer, Berlin, Heidelberg.

Wölfel, M., & Ekenel, H. K. (2005, September). Feature weighted Mahalanobis distance: improved robustness for Gaussian classifiers. In 2005 13th European signal processing conference (pp. 1-4). IEEE.

Wagner, Maggie R, Derek S Lundberg, G Tijana, Susannah G Tringe, Jeffery L Dangl, and Thomas Mitchell-Olds. 2016. “Host Genotype and Age Shape the Leaf and Root Microbiomes of a Wild Perennial Plant.” Nature Communications 7. Nature Publishing Group: 12151. [doi:10.1038/ncomms12151](doi:10.1038/ncomms12151).