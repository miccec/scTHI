#' scTHI_runTsne
#'
#' Runs t-SNE dimensionality reduction on selected features based on Rtsne
#' package.
#' @param scTHIresult scTHI object.
#' @examples
#' library(scTHI.data)
#' data(scExample)
#' result <-  scTHI_score(scExample,
#'                        cellCusterA = colnames(scExample)[1:30],
#'                        cellCusterB = colnames(scExample)[31:100],
#'                        cellCusterAName = "ClusterA",
#'                        cellCusterBName = "ClusterB", filterCutoff = 0,
#'                        pvalueCutoff = 1, nPermu = 100, ncore = 8)
#' result <- scTHI_runTsne(result)
#' @return The same object as scTHI_score with a fifth item tsneData
#'   (data.frame)
#' @export
#'
#' scTHI_runTsne

scTHI_runTsne <- function(scTHIresult) {

  expMat <- scTHIresult$expMat
  
  ## create tsne
  eta <- .1
  filter <- apply(expMat, 1, function(x) {
    sum(quantile(x,
                 probs = c(1 - eta, eta)
    ) * c(1, -1))
  })
  # fivenum(filter)
  # plot(density(log2(filter)))
  cutoff <- density(log2(filter))
  cutoff <- data.frame(x = cutoff$x, y = cutoff$y)
  cutoff <- cutoff$x[which.max(cutoff$y)]
  foldChange <- 2^cutoff
  # sum(filter > foldChange) # 4567
  variableGenes <- names(filter)[filter > foldChange]
  expMat <- expMat[variableGenes, ]
  
  requireNamespace("Rtsne")
  expMatT <- t(expMat)
  tsne_out <- Rtsne::Rtsne(expMatT)
  tsneData <- data.frame(
    x = tsne_out$Y[, 1],
    y = tsne_out$Y[, 2],
    Sample = colnames(expMat),
    stringsAsFactors = FALSE
  )
  rownames(tsneData) <- colnames(expMat)
  scTHIresult <- c(list(tsneData), scTHIresult)
  names(scTHIresult)[1] <- "tsneData"
  return(scTHIresult)
}