#' TME_plot
#'
#' Generates a plot on the t-SNE coordinates, labeling cells by TME
#' classification.
#'
#' @param tsneData X and y coordinates of points in the plot.
#' @param Class Object returned by TME_classification function.
#' @param cexPoint Set the point size.
#' @examples
#' library(scTHI.data)
#' data(scExample)
#' result <-  scTHI_score(scExample,
#'            cellCusterA = colnames(scExample)[1:30],
#'            cellCusterB = colnames(scExample)[31:100],
#'            cellCusterAName = "ClusterA",
#'            cellCusterBName = "ClusterB", filterCutoff = 0,
#'            pvalueCutoff = 1, nPermu = 100, ncore = 8)
#' result <- scTHI_runTsne(result)
#' Class <- TME_classification(scExample)
#' TME_plot(tsneData = result$tsneData, Class)
#' @return  None
#' @export
#'
#' TME_plot

TME_plot <- function(tsneData, Class, cexPoint = 0.8) {
  
  ClassColor <- Class$ClassLegend[Class$Class]
  names(ClassColor) <- names(Class$Class)
  
  TmeColors <- rep("gray80", nrow(tsneData))
  names(TmeColors) <- rownames(tsneData)
  TmeColors[names(ClassColor)] <- ClassColor
  par(oma = c(0, 0, 0, 6))
  plot(
    tsneData[, 1:2],
    pch = 16,
    cex = cexPoint,
    col = TmeColors[rownames(tsneData)],
    main = "",
    xlab = "tSNE 1",
    ylab = "tSNE 2"
  )
  legend(
    par("usr")[2],
    par("usr")[4],
    bty = "n",
    xpd = NA,
    legend = names(Class$ClassLegend),
    col = Class$ClassLegend,
    cex = 0.9,
    pch = 15,
    box.lty = 0
  )
}
