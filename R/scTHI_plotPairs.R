colByValue_ <-
  function(x,
           col,
           range = NA,
           breaks = NA,
           cex.axis = 2,
           las = 1,
           ...) {
    if (is.vector(x)) {
      x <- as.matrix(x)
    }
    if (is.na(range[1])) {
      
    }
    else {
      x[x < range[1]] <- range[1]
      x[x > range[2]] <- range[2]
    }
    if (is.na(breaks[1])) {
      ff <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                length = length(col)
      )
      bg2 <- apply(
        as.matrix(as.numeric(unlist(x))), 1,
        function(x) {
          rank(c(ff, x), ties.method = "min")[length(col) + 1]
        }
      )
      dens <- matrix(bg2, nrow(x), ncol(x))
      result <- matrix(col[dens], nrow = nrow(x), ncol = ncol(x))
      row.names(result) <- row.names(x)
      # image(x = 1:2, y = as.matrix(ff), z = t(ff), col = col,
      #       xaxt = "n", ylab = "", las = las, xlab = "",
      #       xlim = c(1, 4), bty = "n", ...)
      return(result)
    }
    else {
      temp <- cut(as.numeric(unlist(x)),
                  breaks = breaks,
                  include.lowest = TRUE
      )
      if (length(col) != length(levels(temp))) {
        stop("length:col != length: cut result")
      }
      result <- matrix(col[as.numeric(temp)],
                       nrow = nrow(x),
                       ncol = ncol(x)
      )
      row.names(result) <- row.names(x)
      return(result)
    }
  }



getcolors <- function(genesToplot, expMat, 
                      tsneData){
  
  ans <- lapply(genesToplot, function(x){
    tmp_exp <-
      t(expMat[x, colnames(expMat), drop = FALSE])
    mean_exp_samples <- mean(tmp_exp)
    
    ##### range gene x gene
    range_exp <- range(tmp_exp)
    bre <- seq(
      min(tmp_exp[, x]) - 0.01,
      max(tmp_exp[, x]) + 0.01, 0.01
    )
    colori <-
      colByValue_(
        tmp_exp,
        col = colorRampPalette(c(
          "gray80", "coral",
          "red"
        ))(length(bre) - 1),
        breaks = bre,
        cex.axis = 0.8
      )
    colori <- colori[rownames(tsneData), , drop = FALSE]
    tmp <- list(colori, range_exp, mean_exp_samples)
    names(tmp) <- c("colori", "range_exp", "mean_exp_samples")
    tmp})
  
  return(ans)
}


#' scTHI_plotPairs
#'
#' Generates a plot on the t-SNE coordinates to
#' show the expression levels
#' of an interaction pair of interest. Each cell is
#' colored according to the
#' corresponding
#' gene expression value.
#' @param scTHIresult scTHI object.
#' @param cexPoint Set the point size.
#' @param interactionToplot Interaction pair to plot.
#' @examples
#' library(scTHI.data)
#' data(scExample)
#' result <-  scTHI_score(scExample,
#'                  cellCusterA = colnames(scExample)[1:30],
#'                  cellCusterB = colnames(scExample)[31:100],
#'                  cellCusterAName = "ClusterA",
#'                  cellCusterBName = "ClusterB", filterCutoff = 0,
#'                  pvalueCutoff = 1, nPermu = 100, ncore = 8)
#' result <- scTHI_runTsne(result)
#' scTHI_plotPairs(result,interactionToplot = "CXCL12_CD4")

#' @return None
#' @export
#'
#' scTHI_plotPairs
scTHI_plotPairs <-
  function(scTHIresult,
           cexPoint = 0.8,
           interactionToplot) {
   
    tsneData <- scTHIresult$tsneData
    result <- scTHIresult$result
    expMat <- scTHIresult$expMat
    
    genesToplotA <- result[interactionToplot, c("partnerA")]
    genesToplotA <- unlist(strsplit(genesToplotA, ":"))
    genesToplotB <- result[interactionToplot, c("partnerB")]
    genesToplotB <- unlist(strsplit(genesToplotB, ":"))
    
    #################  create list colori
    list_colori_A <- getcolors(genesToplotA, expMat, tsneData)
    names(list_colori_A) <- genesToplotA
    
    list_colori_B <- getcolors(genesToplotB, expMat, tsneData)
    names(list_colori_B) <- genesToplotB
    
    ################# plot
    par(
      mfrow = c(1, length(c(
        genesToplotA, genesToplotB
      ))),
      mar = c(6, 4, 3, 2),
      cex.main = 1.5,
      cex.sub = 1.2,
      col.main = "black",
      col.sub = "gray30",
      font.main = 2,
      font.sub = 3
    )
    
    lapply(names(list_colori_A), function(x) {
      colori <- list_colori_A[[x]]$colori
      range_exp <- round(list_colori_A[[x]]$range_exp, digits = 2)
      mean_exp_samples <-
        round(list_colori_A[[x]]$mean_exp_samples, digits = 2)
      plot(
        tsneData[, 1:2],
        pch = 16,
        cex = cexPoint,
        col = colori,
        main = x,
        xlab = "",
        ylab = "",
        sub = paste0(
          "Range: ",
          range_exp[1],
          " to ",
          range_exp[2],
          ", MeanExp: ",
          mean_exp_samples
        )
      )
    })
    
    lapply(names(list_colori_B), function(x) {
      colori <- list_colori_B[[x]]$colori
      range_exp <- round(list_colori_B[[x]]$range_exp, digits = 2)
      mean_exp_samples <-
        round(list_colori_B[[x]]$mean_exp_samples, digits = 2)
      plot(
        tsneData[, 1:2],
        pch = 16,
        cex = cexPoint,
        col = colori,
        main = x,
        xlab = "",
        ylab = "",
        sub = paste0(
          "Range: ",
          range_exp[1],
          " to ",
          range_exp[2],
          ", MeanExp: ",
          mean_exp_samples
        )
      )
    })
  
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  }