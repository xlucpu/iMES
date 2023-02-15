# iMES

<!-- badges: start -->
<!-- badges: end -->

## Introduction

This package provides a function to calculate an index of methylation-based epigenetic silencing (iMES) using binary DNA methylation status for patients with clear cell renal cell carcinoma.

## Citation

As iMES paper has not been published, if you use iMES R package in published research, for now please cite:

  - Lu X, Vano Y, Helleux A, Su X, Lindner V, Davidson G, Mouawad R, Spano JP, Roupret M, Elaidi R, Comp\xE9rat E, Verkarre V, Sun C, Chevreau C, Bennamoun M, Lang H, Tricard T, Cheng W, Xu L, Davidson I, Yan F, Fridman WH, Sautes-Fridman C, Oudard S, Malouf GG. An enhancer demethylator phenotype converged to immune dysfunction and resistance to immune checkpoint inhibitors in clear-cell renal cell carcinomas. Clin Cancer Res. 2022 Nov 14:CCR-22-2133. doi: 10.1158/1078-0432.CCR-22-2133.

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
methMat <- read.table("DNA methylation beta matrix.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
iMES <- iMES(bmat     = methMat, # a DNA methylation beta matrix with continuous values as input
             methcut  = 0.2, # cut continuous methylation matrix to binary methylation status
             samples  = colnames(methMat)[1:30], # extract the first 30 samples to calculate iMES
             quantile = 3) # dichotomize samples into iMES-high and iMES-low based on a general tertile cutoff
```

