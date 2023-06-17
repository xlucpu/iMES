# iMES
<!-- badges: start -->
<!-- badges: end -->

## Table of Contents
1. [Introduction](#introduction)
2. [Citation](#citation)
3. [Installation](#installation)
4. [Example](#example)

## Introduction <a name="introduction"></a>
This package is designed to compute an Index of Methylation-based Epigenetic Silencing (iMES) using binary DNA methylation data in patients with clear cell renal cell carcinoma. Furthermore, it classifies patients into distinct regulon phenotypes based on transcriptomic expression data.

The primary function, `iMES()`, returns a DataFrame where each row corresponds to a sample name and comprises three columns:
1. 'iMES': raw iMES score
2. 'iMES.mm': min-max normalized iMES score (scaled to a range of 0-10) multiplied by 10
3. 'iMES.group': dichotomized iMES group

The secondary function, `predRegulon()`, assesses the regulon activity status for each sample. It then categorizes samples into either a suppressed or activated regulon phenotype, based on the status count of each regulon within each original group.

Patients with high iMES scores or those categorized in the iMES-high group are potentially at an increased risk of immune evasion and resistance to immune checkpoint inhibitors. The regulon phenotype classified as suppressed is considered to correspond to the iMES-high group, while an activated regulon phenotype is analogous to the iMES-low group.

## Citation <a name="citation"></a>
As iMES paper has not been published, if you use iMES R package in published research, for now please cite:

  - Lu X, Vano Y, Helleux A, Su X, Lindner V, Davidson G, Mouawad R, Spano JP, Roupret M, Elaidi R, Comp&eacute;rat E, Verkarre V, Sun C, Chevreau C, Bennamoun M, Lang H, Tricard T, Cheng W, Xu L, Davidson I, Yan F, Fridman WH, Sautes-Fridman C, Oudard S, Malouf GG. An enhancer demethylator phenotype converged to immune dysfunction and resistance to immune checkpoint inhibitors in clear-cell renal cell carcinomas. Clin Cancer Res. 2022 Nov 14:CCR-22-2133. doi: 10.1158/1078-0432.CCR-22-2133.

## Installation <a name="installation"></a>
You can install the development version of iMES from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xlucpu/iMES")