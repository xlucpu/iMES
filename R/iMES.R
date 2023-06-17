#' @name iMES
#' @title Index of methylation-based epigenetic silencing
#' @description This function calculates an index of methylation-based epigenetic silencing (iMES) using binary DNA methylation status for patients with clear cell renal cell carcinoma.
#' @author Xiaofan Lu
#'
#' @param bmat A numeric DNA methylation beta matrix with row features (probes) and sample columns and continuous values as input.
#' @param methcut A numeric value to indicate the methylation cutoff and assign each probe to be either methylated or unmethylated; 0.2 by default.
#' @param samples A string value to indicate the samples that will be used to calculate iMES; all samples will be used by default.
#' @param quantile A numeric value to indicate quantile base to dichotomize samples into iMES-high and iMES-low; 3 (tertile) by default.
#'
#' @return A DataFrame with rownames of samples and three columns: iMES (raw iMES score), iMES.mm (minmax normalized iMES score * 10; range from 0-10), iMES.group (dichotomized iMES group)
#' @export
#' @importFrom lsr quantileCut
#'
iMES <- function(bmat     = NULL,
                 methcut  = 0.2,
                 samples  = NULL,
                 quantile = 3) {

  # customized function for min-max normalization
  range01 <- function(x){(x - min(x)) / (max(x) - min(x)) * 10}

  # preprocess data
  bmat <- as.data.frame(na.omit(bmat))

  # check if all probes can be matched
  if(all(is.element(adaLASSO.coeff$probe, rownames(bmat)))) {
    bmat <- bmat[adaLASSO.coeff$probe, , drop = F]
  } else {
    missPb <- setdiff(adaLASSO.coeff$probe, rownames(bmat))
    if(length(missPb) == 1) {
      stop("missing the following probe: ", missPb)
    } else {
      stop("missing the following probes: ", paste(missPb, collapse = ", "))
    }
  }

  # transfer continous DNA methylation data to binary status using methylation cutoff
  bmat[bmat >  methcut] <- 1
  bmat[bmat <= methcut] <- 0

  # calculate iMES
  if(is.null(samples)) {
    samples  <- colnames(bmat)
    iMES.raw <- apply(t(bmat[adaLASSO.coeff$probe, ,drop = FALSE]), 1, function(x) {x %*% adaLASSO.coeff$coeff})
  } else {
    iMES.raw <- apply(t(bmat[adaLASSO.coeff$probe, samples, drop = FALSE]), 1, function(x) {x %*% adaLASSO.coeff$coeff})
  }

  # if normalize
  iMES.mm <- range01(iMES.raw)

  # generate data frame
  iMES <- data.frame(iMES = as.numeric(iMES.raw),
                     iMES.minmax = as.numeric(iMES.mm),
                     row.names = samples)

  # dichotomize samples
  cut <- quantileCut(iMES$iMES,quantile)
  iMES.group <- ifelse(as.character(cut) == levels(cut)[quantile],"iMES-high","iMES-low")
  iMES$iMES.group <- iMES.group

  # return
  return(iMES)
}
