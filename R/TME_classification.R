#' TME_classification
#'
#' The function allows the user to classify non-tumor cells in tumor
#' microenvironment. It implements the Mann-Whitney-Wilcoxon Gene
#' Set Test
#' (MWW-GST) algorithm and tests for each cell the enrichment of
#' a collection
#' of signatures of different cell types.
#' @param expMat Gene expression matrix where rows are genes
#'   presented with
#'   Hugo Symbols and columns are cells. Gene expression values
#'   should be normalized counts.
#' @param minLenGeneSet Minimum gene set length
#' @param pvalFilter Logical, if TRUE results will be filtered
#' for p-Value.
#'   Defoult is FALSE.
#' @param alternative a character string specifying the alternative
#' hypothesis of wilcoxon test, must be one of "two.sided" (default),
#' "greater" or "less".
#' @param fdrFilter Logical, if TRUE results will be filtered for FDR.
#' @param pvalCutoff Numeric p-Value (or FDR) threshold. Gene set with
#'   p-Value (or FDR) greater than pvalCutoff will be discarded
#'   (default is 0.01).
#' @param nesCutoff Numeric threshold. Gene set with NES greater than
#'   nesCutoff will be discarded (default is 0.58)
#' @param nNES Default is 0.58, so each cell is classified with
#'  a specific
#'   phenotype based on the first significant enriched gene set.
#' @examples
#' library(scTHI.data)
#' data(scExample)
#' Class <- TME_classification(scExample)
#' @return A list with two items: Class (character) and ClassLegend
#'   (character)
#' @export
#'
#' TME_classification

TME_classification <- function(expMat,
                               minLenGeneSet = 10,
                               alternative = "two.sided",
                               pvalFilter = FALSE,
                               fdrFilter = TRUE,
                               pvalCutoff = 0.01,
                               nesCutoff = 0.58,
                               nNES = 1) {
  
  zerogenes <- rowSums(expMat == 0)
  expMat <- expMat[zerogenes != ncol(expMat), ]
  means <- rowMeans(expMat)
  sds <- apply(expMat, 1, sd)
  E <- apply(expMat, 2, function(x) {
    (x - means) / sds
  })
  
  E_ <- apply(E, 2, rank)
  tmp <- lapply(signatures, function(x) {
    sum(x %in% rownames(E))
  })
  signatures_ <- signatures[tmp > minLenGeneSet]
  
  ans <- vector("list", length(signatures_))
  names(ans) <- names(signatures_)
  ans <- lapply(signatures_, function(S) {
    geneSet <- S
    gs <- geneSet[which(geneSet %in% rownames(E))]
    outside_gs <- setdiff(rownames(E), gs)
    nx <- length(gs)
    ny <- length(outside_gs)
    Tstat <- colSums(E_[outside_gs, ])
    Ustat <- nx * ny + ny * (ny + 1) / 2 - Tstat
    mu <- nx * ny / 2
    sigma <- sqrt(mu * (nx + ny + 1) / 6)
    zValue <- Ustat - mu
    correction <- vapply(zValue, function(x) {
      switch(
        alternative,
        two.sided = sign(x) * 0.5,
        greater = 0.5,
        less = -0.5
      )
    },numeric(1))
    zValue <- (zValue - correction) / sigma
    pValue <- vapply(zValue, function(x) {
      switch(
        alternative,
        less = 1 - pnorm(x),
        greater = pnorm(-x),
        two.sided = 2 * pnorm(-abs(x))
      )
    },numeric(1))
    nes <- Ustat / nx / ny
    pu <- nes / (1 - nes)
    log.pu <- log2(pu)
    names(Ustat) <- colnames(E_)
    names(pValue) <- colnames(E_)
    names(nes) <- colnames(E_)
    names(pu) <- colnames(E_)
    names(log.pu) <- colnames(E_)
    #ans[[i]] <- 
    list(
      statistic = Ustat,
      p.value = pValue,
      nes = nes,
      pu = pu,
      log.pu = log.pu
    )})
  
  
  NES <- vapply(ans, function(x) {
    cbind(x$log.pu)
  },numeric(ncol(expMat)))
  NES <- t(NES)
  colnames(NES) <- colnames(expMat)
  pValue <- vapply(ans, function(x) {
    cbind(x$p.value)
  },numeric(ncol(expMat)))
  pValue <- t(pValue)
  colnames(pValue) <- colnames(expMat)
  
  FDR <- apply(pValue, 2, function(x) {
    p.adjust(x, method = "fdr")
  })
  if (pvalFilter == TRUE) {
    NES[pValue >= pvalCutoff] <- 0
  }
  if (fdrFilter == TRUE) {
    NES[FDR >= pvalCutoff] <- 0
  }
  NES[NES < nesCutoff] <- 0
  if (nNES == 1) {
    Class <- apply(NES, 2, function(x) {
      names(which.max(x))
    })
  }
  if (nNES != 1) {
    Class <- apply(NES, 2, function(x) {
      names(sort(x, decreasing = TRUE))[nNES]
    })
  }
  Class[colSums(NES != 0) < nNES] <- "nc"
  phenotype <- signaturesColors[Class, ]
  rownames(phenotype) <- names(Class)
  Class <- phenotype$ALLPhenotypeFinal
  names(Class) <- rownames(phenotype)
  ClassLegend <- phenotype$Color
  names(ClassLegend) <- phenotype$ALLPhenotypeFinal
  ClassLegend <- ClassLegend[!duplicated(ClassLegend)]
  #print(sort(table(Class), decreasing = TRUE))
  Classification <- list(Class, ClassLegend)
  names(Classification) <- c("Class", "ClassLegend")
  return(Classification)
}