#' AdaLASSO coefficient
#'
#' A data frame storing feature coefficient for model-selected probes
#'
#' @format A data frame with 58 rows (probes relevant to silenced genes) and their corresponding adaLASSO coefficient
#' \describe{contains coefficient to calculate iMES.}
"adaLASSO.coeff"

#' TNI object
#'
#' An R object derived from Transcriptional Network Inference
#'
#' @format An object of class TNI
#' \describe{contains TNI an object of class RTN and can be used for external prediction of regulon activity.}
"rtni_kirc"

#' mRNA list
#'
#' An vector including mRNAs
#'
#' @format An vector including mRNAs
#' \describe{contains 19,620 mRNAs.}
"Mids"

#' LASSO gene features
#'
#' An vector of genes that constitutes iMES
#'
#' @format An vector of genes that constitutes iMES
#' \describe{contains 55 genes.}
"lasso_fea_gene"
