#' scTHI_plotCluster
#'
#' Graphs the output of scTHI_runTsne, labeling cells by clusters.
#' @param scTHIresult scTHI object.
#' @param cexPoint Set the point size.
#' @param legendPos Character string to custom the legend position.
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
#' scTHI_plotCluster(result)
#' @return None
#' @export
#'
#' scTHI_plotCluster

scTHI_plotCluster <- function(scTHIresult,
                              cexPoint = 0.8,
                              legendPos = c(
                                "topleft", "topright",
                                "bottomright", "bottomleft"
                              )) {
  
  tsneData <- scTHIresult$tsneData
  legendPos <- match.arg(legendPos)
  
  nCluster <- length(scTHIresult) - 3
  ggplot_colors <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  Colors <- ggplot_colors(n = nCluster)
  names(Colors) <- names(scTHIresult)[4:length(scTHIresult)]
  
  ClusterColors <- rep("gray80", nrow(tsneData))
  names(ClusterColors) <- rownames(tsneData)
  
  ClusterColors[scTHIresult[[names(Colors)[1]]]] <- Colors[1]
  ClusterColors[scTHIresult[[names(Colors)[2]]]] <- Colors[2]
  
  plot(
    tsneData[, 1:2],
    pch = 16,
    cex = cexPoint,
    col = ClusterColors[rownames(tsneData)],
    main = "",
    xlab = "tSNE 1",
    ylab = "tSNE 2"
  )
  legend(
    legendPos,
    inset = .02,
    legend = names(Colors),
    col = Colors,
    cex = 1,
    pch = 15,
    box.lty = 0,
    bg = "transparent"
  )
}