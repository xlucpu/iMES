# iMES

<!-- badges: start -->
<!-- badges: end -->

## Introduction

This package provides a function for calculating an index of methylation-based epigenetic silencing (iMES) using binary DNA methylation status for patients with clear cell renal cell carcinoma. Additionally, it enables the classification of patients into specific regulon phenotypes using transcriptomic expression data. The main function, iMES(), returns a DataFrame with the names of the samples as rows and three columns: iMES (raw iMES score), iMES.mm (min-max normalized iMES score multiplied by 10, with a range of 0-10), and iMES.group (a dichotomized iMES group). The second function, predRegulon, calculates regulon activity status for each sample and classify samples into suppressed or activated regulon phenotype by counting the status of each regulon for each original group. Patients with high iMES scores or those categorized in the iMES-high group may have a higher likelihood of developing immune invasion and resistance to immune checkpoint inhibitors, and regulon phenotype that was classified as suppressed was considered to mirror iMES-high group, while activated regulon phenotype represented iMES-low group.
## Citation

As iMES paper has not been published, if you use iMES R package in published research, for now please cite:

  - Lu X, Vano Y, Helleux A, Su X, Lindner V, Davidson G, Mouawad R, Spano JP, Roupret M, Elaidi R, Comp&eacute;rat E, Verkarre V, Sun C, Chevreau C, Bennamoun M, Lang H, Tricard T, Cheng W, Xu L, Davidson I, Yan F, Fridman WH, Sautes-Fridman C, Oudard S, Malouf GG. An enhancer demethylator phenotype converged to immune dysfunction and resistance to immune checkpoint inhibitors in clear-cell renal cell carcinomas. Clin Cancer Res. 2022 Nov 14:CCR-22-2133. doi: 10.1158/1078-0432.CCR-22-2133.

## Installation

You can install the development version of iMES from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xlucpu/iMES")
```

## Example
``` r
## basic example code (not run)
library(iMES)
methMat <- read.delim("DNA methylation beta matrix.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
iMES <- iMES(bmat     = methMat, # a DNA methylation beta matrix with continuous values as input
             methcut  = 0.2, # cut continuous methylation matrix to binary methylation status
             samples  = colnames(methMat)[1:30], # extract the first 30 samples to calculate iMES
             quantile = 3) # dichotomize samples into iMES-high and iMES-low based on a general tertile cutoff
			 
exprMat <- 	read.delim("Transcriptomic expression matrix.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)	
regulon <- predRegulon(emat     = exprMat,
                       samples  = colnames(exprMat)[1:30], # extract the first 30 samples to calculate regulon activity
                       seed     = 20000112, # seed to reproduce K-mode clustering (k = 2)
                       fig.path = getwd(), # path to save regulon activity heatmap
                       fig.name = "heatmap of regulon activity") # name of the regulon activity heatmap
```