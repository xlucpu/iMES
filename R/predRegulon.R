#' @name predRegulon
#' @title Regulon phenotypes based on genes with epigenetic silencing
#' @description This function inferrs regulon activity of genes that were epigenetically silenced using transcriptomic expression data for patients with clear cell renal cell carcinoma.
#' @author Xiaofan Lu
#'
#' @param emat A numeric transcritomic gene expression with row features (genes) and sample columns and continuous values as input with proper normalisation (e.g., TPM, FPKM, normalized count, or microarray signals).
#' @param samples A string value to indicate the samples that will be used to calculate regulon activity; all samples will be used by default.
#' @param seed A numeric string to indicate seed for K-mode clustering for reproducibility
#' @param fig.path A string value to indicate the output path for storing the regulon activity heatmap.
#' @param fig.name A string value to indicate the name of the regulon activity heatmap.
#' @return A DataFrame with rownames of regulons and colnames of samples with input of regulon activity status, and a predictive regulon phenotype based on K-modes clustering (k=2) with a heatmap.
#' @export
#' @importFrom RTN tni.replace.samples tni.gsea2 tni.get
#' @importFrom klaR kmodes
#' @importFrom ClassDiscovery distanceMatrix
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom grid unit
#'
predRegulon <- function(emat     = NULL,
                        samples  = NULL,
                        seed     = 20000112,
                        fig.path = getwd(),
                        fig.name = "heatmap of regulon activity") {

  # customized function
  quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
      sink(tempfile())
      on.exit(sink())
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
  }

  get_status <- function(df, group) {
    df_group <- df[, group]
    c(sum(df_group == 1), sum(df_group == -1))
  }

  # emat <- read.delim("expr.kirc.txt", sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  # load("rtni_kirc.RData")
  # load("Mids.RData")
  # load("lasso_fea_gene.RData")

  # check scale of the expression data
  if(max(emat) < 25 | (max(emat) >= 25 & min(emat) < 0)) {
    message("--expression profile seems to have veen standardised (z-score or log transformation), no more action will be performed.")
    gset <- emat
  }
  if(max(emat) >= 25 & min(emat) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(emat + 1)
  }

  if(is.null(samples)) {
    samples <- colnames(gset)
  }

  # calculate regulon activity in external dataset
  rtni_myDT <- quiet(tni.replace.samples(rtni_kirc, as.matrix(gset[intersect(Mids,rownames(gset)),samples])))
  rtnigsea_myDT <- quiet(tni.gsea2(rtni_myDT, regulatoryElements = intersect(rownames(gset),lasso_fea_gene)))
  myDT_regact <- tni.get(rtnigsea_myDT, what = "regulonActivity")

  # extract regulon activity
  indata <- as.data.frame(t(myDT_regact$status))
  indata <- indata[,samples]

  # perform K-mode clustering
  set.seed(seed)
  kmodes_result <- kmodes(t(indata), 2)
  hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")

  # generate annotation
  annCol <- data.frame(row.names = samples,
                       Cluster = paste0("C",kmodes_result$cluster))
  annColors <- list("Cluster" = c("C1" = "#7CDDD2","C2" = "#B21558"),
                    "Regulon" = c("Suppressed" = "#7CDDD2","Activated" = "#B21558"))
  col_fun <- colorRamp2(c(-1, 0, 1), c("#7CDDD2", "white", "#B21558"))

  # determine regulon-activated and regulon-suppressed subtypes
  res <- paste0("C",kmodes_result$cluster); names(res) <- samples
  C1 <- names(res)[res == "C1"]
  C2 <- names(res)[res == "C2"]

  # Loop through the genes
  df <- indata
  for (gene in rownames(df)) {

    # Get status of the gene in C1 and C2
    status_C1 <- get_status(df[gene, ], C1)
    status_C2 <- get_status(df[gene, ], C2)

  }

  # Classify the groups
  C1_active <- sum(df[, C1] == 1)
  C2_active <- sum(df[, C2] == 1)

  if (C1_active > C2_active) {
    print("C1 is globally activated")
    annCol$Regulon <- ifelse(annCol$Cluster == "C1","Activated","Suppressed")
    annCol$Regulon <- factor(annCol$Regulon, levels = c("Suppressed","Activated"))
    samorder <- rownames(annCol[order(annCol$Regulon),,drop = F])

    top_anno <- HeatmapAnnotation(df = annCol[samorder,"Regulon",drop = F],
                                  col = annColors,
                                  height = unit(10,"mm"),
                                  border = T)

    hm <- Heatmap(as.matrix(indata[,samorder]),
                  col = col_fun,
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  cluster_rows = hcg,
                  cluster_columns = F,
                  heatmap_height = unit(300, "cm")/nrow(indata),
                  heatmap_width = unit(8, "cm"),
                  top_annotation = top_anno,
                  heatmap_legend_param = list(
                    at = c(-1, 0, 1),
                    labels = c("Inactive", "Wild", "Active"),
                    title = "Activity",
                    color_bar = "discrete"))
    pdf(file = file.path(fig.path,paste0(fig.name, ".pdf")), width = 6,height = 6)
    draw(hm, heatmap_legend_side = "right",annotation_legend_side = "right")
    invisible(dev.off())

  } else if (C1_active < C2_active) {
    print("C2 is globally activated")
    annCol$Regulon <- ifelse(annCol$Cluster == "C2","Activated","Suppressed")
    annCol$Regulon <- factor(annCol$Regulon, levels = c("Suppressed","Activated"))
    samorder <- rownames(annCol[order(annCol$Regulon),,drop = F])

    top_anno <- HeatmapAnnotation(df = annCol[samorder,"Regulon",drop = F],
                                  col = annColors,
                                  height = unit(10,"mm"),
                                  border = T)

    hm <- Heatmap(as.matrix(indata[,samorder]),
                  col = col_fun,
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  cluster_rows = hcg,
                  cluster_columns = F,
                  heatmap_height = unit(300, "cm")/nrow(indata),
                  heatmap_width = unit(8, "cm"),
                  top_annotation = top_anno,
                  heatmap_legend_param = list(
                    at = c(-1, 0, 1),
                    labels = c("Inactive", "Wild", "Active"),
                    title = "Activity",
                    color_bar = "discrete"))
    pdf(file = file.path(fig.path,paste0(fig.name, ".pdf")), width = 6,height = 6)
    draw(hm, heatmap_legend_side = "right",annotation_legend_side = "right")
    invisible(dev.off())

  } else {
    print("C1 and C2 have an equal number of global activations")
    annCol$Cluster <- factor(annCol$Cluster, levels = c("C1","C2"))
    samorder <- rownames(annCol[order(annCol$Cluster),,drop = F])

    top_anno <- HeatmapAnnotation(df = annCol[samorder,"Cluster",drop = F],
                                  col = annColors,
                                  height = unit(10,"mm"),
                                  border = T)

    hm <- Heatmap(as.matrix(indata[,samorder]),
                  col = col_fun,
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  cluster_rows = hcg,
                  cluster_columns = F,
                  heatmap_height = unit(300, "cm")/nrow(indata),
                  heatmap_width = unit(8, "cm"),
                  top_annotation = top_anno,
                  heatmap_legend_param = list(
                    at = c(-1, 0, 1),
                    labels = c("Inactive", "Wild", "Active"),
                    title = "Activity",
                    color_bar = "discrete"))
    pdf(file = file.path(fig.path,paste0(fig.name, ".pdf")), width = 6,height = 6)
    draw(hm, heatmap_legend_side = "right",annotation_legend_side = "right")
    invisible(dev.off())
  }

  return(list(Status = as.data.frame(t(myDT_regact$status)),
              Phenotype = annCol))
}
