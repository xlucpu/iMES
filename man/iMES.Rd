% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iMES.R
\name{iMES}
\alias{iMES}
\title{Index of methylation-based epigenetic silencing}
\usage{
iMES(bmat = NULL, methcut = 0.2, samples = NULL, quantile = 3)
}
\arguments{
\item{bmat}{A numeric DNA methylation beta matrix with row features (probes) and sample columns and continuous values as input.}

\item{methcut}{A numeric value to indicate the methylation cutoff and assign each probe to be either methylated or unmethylated; 0.2 by default.}

\item{samples}{A string value to indicate the samples that will be used to calculate iMES; all samples will be used by default.}

\item{quantile}{A numeric value to indicate quantile base to dichotomize samples into iMES-high and iMES-low; 3 (tertile) by default.}
}
\value{
A DataFrame with rownames of samples and three columns: iMES (raw iMES score), iMES.mm (minmax normalized iMES score * 10; range from 0-10), iMES.group (dichotomized iMES group)
}
\description{
This function calculates an index of methylation-based epigenetic silencing (iMES) using binary DNA methylation status for patients with clear cell renal cell carcinoma.
}
\author{
Xiaofan Lu
}
