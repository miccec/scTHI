#' scTHI_plotResult
#'
#' Creates barplots of scTHI_score results.
#' @param scTHIresult scTHI object.
#' @param cexNames Size of names in barplot.
#' @param plotType Type of plot to be generated. Default is
#' "score", can be  also "pair". The "score" option will generate
#' a barplot for each resulted interaction pair, representing
#' the calculated interaction score
#'   and the related p-Value.The "pair" option will generate
#'   two barplot for each resulted interaction pair, representing
#'   the percentage of cells of each cluster expressing partnerA
#'   and partnerB gene, respectively.
#' @param nRes Number of pairs to plot (all if NULL).
#' @examples
#' data(scExample)
#' result <-  scTHI_score(scExample,
#'                        cellCusterA = colnames(scExample)[1:30],
#'                        cellCusterB = colnames(scExample)[31:100],
#'                        cellCusterAName = "ClusterA",
#'                        cellCusterBName = "ClusterB", filterCutoff = 0,
#'                        pvalueCutoff = 1, nPermu = 100, ncore = 8)
#'
#' scTHI_plotResult(result, plotType = "score")
#' scTHI_plotResult(result, plotType = "pair")
#' @return None
#' @export
#'
#' scTHI_plotResult

scTHI_plotResult <- function(scTHIresult,
                             cexNames = 0.8,
                             plotType = c("score", "pair"),
                             nRes = NULL) {
  
  result <- scTHIresult$result
  if (!is.null(nRes)) {
    result <- result[1:nRes, ]
  }
  
  
  if (plotType == "score") {
    tmp <- as.matrix(result[, c("SCORE", "pValue")])
    pvalues <- apply(tmp, 1, function(x) {
      if (x["pValue"] > 0.05) "ns" else if (x["pValue"] <= 0.0001){ 
        "****"
      }else if (x["pValue"] <= 0.001) {
        "***"
      }else if (x["pValue"] <= 0.01) {
        "**"
      }else if (x["pValue"] <= 0.05) {
        "*"
      }
      
    })
    
    par(mar = c(5, 10, 4, 2) + 0.1)
    barplot(
      t(tmp[, "SCORE"]),
      beside = TRUE,
      xlim = c(0, 1.1),
      xlab = "scTHI Score",
      col = c("lightseagreen"),
      cex.names = cexNames,
      horiz = TRUE,
      las = 2
    )
    text(
      x = 1.05,
      y = seq(1.5, by = 2, length.out = nrow(tmp)),
      labels = pvalues,
      cex = 1.2
    )
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
  
  if (plotType == "pair") {
    tmp <- as.matrix(result[, grep("rnk", colnames(result))])
    par(mar = c(5, 10, 4, 2) + 0.1, oma = c(0, 0, 0, 6))
    barplot(
      t(tmp),
      beside = TRUE,
      xlim = c(0, 1.1),
      xlab = "% of Cells",
      col = c("#F8766D", "#00BFC4"),
      cex.names = cexNames,
      horiz = TRUE,
      las = 2
    )
    legend(
      par("usr")[2],
      par("usr")[4],
      bty = "n",
      xpd = NA,
      legend = c("ClusterA", "ClusterB"),
      col = c("#F8766D", "#00BFC4"),
      cex = 1,
      pch = 15,
      box.lty = 0,
      bg = "transparent"
    )
    par(mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
  }
}