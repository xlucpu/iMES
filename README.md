# iMES

<!-- badges: start -->
<!-- badges: end -->

## Introduction

This package provides a function to calculates an index of methylation-based epigenetic silencing (iMES) using binary DNA methylation status for patients with clear cell renal cell carcinoma.

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
methMat <- read.table("DNA methylation beta matrix.txt",,sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
iMES <- iMES(bmat     = methMat, # a DNA methylation beta matrix with continuous values as input
             methcut  = 0.2, # cut continuous methylation matrix to binary methylation status
             samples  = colnames(methMat)[1:30], # extract the first 30 samples to calculate iMES
             quantile = 3) # dichotomize samples into iMES-high and iMES-low based on a general tertile cutoff
```

