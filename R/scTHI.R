
scTHI.score <- function(expMat, cellCusterA, cellCusterB, cellCusterAName, cellCusterBName, topRank = 10, fileNameBase = "scTHI", filterCutoff = 0.50,
                        PValue = TRUE, pvalueCutoff = 0.05, nPermu = 1000, ncore = 8){

  #' scTHI.score
  #'
  #' This function allows the user to compute a score for a set of ligand-receptor pairs, from a single cell gene expression matrix,
  #' and detect specific Tumor-Host interactions. You must specify at least two clusters of cells (for example tumor cells and immune cells).
  #' @param expMat ScRNA-seq gene expression matrix where rows are genes presented with Hugo Symbols and columns are cells. Gene expression values should be counts or normalized counts.
  #' @param clusterA Vector of columns of expMat that belong to the first cluster.
  #' @param clusterB Vector of columns of expMat that belong to the second cluster.
  #' @param cellCusterAName A character string labeling the clusterA.
  #' @param cellCusterBName A character string labeling the clusterB.
  #' @param topRank Filter threshold. Set to 10 (default) means that each gene of the interaction pair will be considered as expressed in a cell if it's in the top rank 10%.
  #' @param fileNameBase Project name.
  #' @param filterCutoff Score threshold (default is 0.50). For each interaction pair, if the score calculated (for the partnerA or partnerB) will be less than filterCutoff the interaction pair will be discarded.
  #' @param foldChange Logical, default is TRUE. This argument allows to calculate a fold change for each gene of the interaction pair (PartenerA and PartnerB) with respect to a third cluster of cells, specified in cellsToCompare.
  #' @param useData Logical, default is FALSE. Set to TRUE allows the user to compute fold change and pvalue later, providing as input a score matrix previously calculated.
  #' @param dataToUse Score matrix to be provided as input if dataToUse is TRUE.
  #' @param cellsToCompare Vector of columns of expMatToCompare that belong to the third cluster.
  #' @param cellsToCompareNAME A character string labeling the third cluster.
  #' @param expMatToCompare ScRNA-seq gene expression matrix of the third cluster cells. It may be the same as expMat.
  #' @param PValue Logical, set to TRUE (default) compute statistical iterations. If p.value < 0.05, the value will be returned.
  #' @param nPermu Number of iterations to perform (default is 1000).
  #' @param ncore Number of processors to use.
  #' @keywords interaction
  #' @export
  #' @examples
  #'
  #' scTHI
    ####################### check rownames expMat
  #load("/storage/gluster/vol1/SHARED/HOMEFOLDERS/caruso/scProject/NYnontumor/Phone/PhoneCustom/package_data/interactionTable/interaction_tableComplete_1144.RData") #### da rimuovere x il pacchetto
  tmp_check <- interaction_table[, c("partnerA1" ,"partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")]
  tmp_check <- unique(tmp_check[!is.na(tmp_check)])
  if(sum(tmp_check %in% rownames(expMat)) != length(tmp_check)){
    tmp_genes <- setdiff(tmp_check, rownames(expMat))
    tmp <- which(interaction_table$partnerA1 %in% tmp_genes | interaction_table$partnerA2 %in% tmp_genes | interaction_table$partnerA3 %in% tmp_genes | interaction_table$partnerB1 %in% tmp_genes | interaction_table$partnerB2 %in% tmp_genes | interaction_table$partnerB3 %in% tmp_genes)
    interaction_table <- interaction_table[-tmp, ]
    message("Warning: Not all interaction genes are in expMat")
  }
  if(nrow(interaction_table) == 0) stop("ERROR: No interaction genes to test")

  ######################### compute score ##################
  message(paste("Computing score for", nrow(interaction_table), "interaction pairs..."))

  ddataA <- expMat[, cellCusterA]
  n <- round(nrow(ddataA) * topRank/100)
  ddataA_ <- apply(ddataA, 2, function(x) rank(-x))
  ddataA_ <- apply(ddataA_, 2, function(x) x <= n)
  rnkA <- rep(0, nrow(interaction_table))
  expValueA <- rep(0, nrow(interaction_table))
  names(rnkA) <- names(expValueA) <- rownames(interaction_table)
  for(i in 1:nrow(interaction_table)){
    ggenes <- interaction_table[i, c("partnerA1", "partnerA2", "partnerA3")]
    ggenes <- as.vector(ggenes[!is.na(ggenes)])
    tmp_genes <- ddataA_[ggenes, , drop = F]
    tmp <- colSums(tmp_genes)
    rnkA[i] <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
    expValueA[i] <- paste(round(rowMeans(ddataA[ggenes, cellCusterA, drop = F]), digits = 2), collapse = " # ")
  }
  resultA <- data.frame(rnkPartnerA = rnkA, expValueA = expValueA, stringsAsFactors = F)
  colnames(resultA) <- paste0(colnames(resultA), "_", cellCusterAName)
  print(paste("Computed", i , "ranked values for partner A"))

  ddataB <- expMat[, cellCusterB]
  n <- round(nrow(ddataB) * topRank/100)
  ddataB_ <- apply(ddataB, 2, function(x) rank(-x))
  ddataB_ <- apply(ddataB_, 2, function(x) x <= n)
  rnkB <- rep(0, nrow(interaction_table))
  expValueB <- rep(0, nrow(interaction_table))
  names(rnkB) <- names(expValueB) <- rownames(interaction_table)
  for(i in 1:nrow(interaction_table)){
    ggenes <- interaction_table[i, c("partnerB1", "partnerB2", "partnerB3")]
    ggenes <- as.vector(ggenes[!is.na(ggenes)])
    tmp_genes <- ddataB_[ggenes, , drop = F]
    tmp <- colSums(tmp_genes)
    rnkB[i] <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
    expValueB[i] <- paste(round(rowMeans(ddataB[ggenes, cellCusterB, drop = F]), digits = 2), collapse = " # ")
  }
  resultB <- data.frame(rnkPartnerB = rnkB, expValueB = expValueB, stringsAsFactors = F)
  colnames(resultB) <- paste0(colnames(resultB), "_", cellCusterBName)
  print(paste("Computed", i , "ranked values for partner B"))

  SCORE <- rowMeans(cbind(rnkA, rnkB))
  result <- data.frame(interaction_table, resultA, resultB, SCORE)

  #### remove low rank-score interaction
  message("Removing low rank-score interactions...")
  interestColumn <- c(paste0("rnkPartnerA_", cellCusterAName), paste0("rnkPartnerB_", cellCusterBName))
  tmp <- rowSums(result[, interestColumn] < filterCutoff)
  result[tmp != 0, "SCORE"] <- NA
  result <- result[!is.na(result$SCORE), ]

  #### prepare final result
  result <- result[order(result$SCORE, decreasing = T), ]
  result <- data.frame(interationPair = rownames(result), result, stringsAsFactors = F)
  columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
  result <- result[, setdiff(colnames(result), columnToremove)]

  if(nrow(result) == 0) stop("No interaction pair exceed the score filterCutoff")
  scTHIresult <- list(result, expMat, cellCusterA, cellCusterB)
  names(scTHIresult) <- c("result", "expMat", cellCusterAName, cellCusterBName)
  save(scTHIresult, file = paste0(fileNameBase, "_", cellCusterAName, "VS", cellCusterBName, ".RData"))


  ######################### compute permutations ##########
  if(PValue == TRUE){
    message("Computing permutation....")

    #load("/storage/gluster/vol1/SHARED/HOMEFOLDERS/caruso/scProject/NYnontumor/Phone/PhoneCustom/package_data/interactionTable/interaction_tableComplete_1144.RData") #### da rimuovere per il pacchetto
    interaction_table <- interaction_table[rownames(result), ]
    columnToadd <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- data.frame(interaction_table[, columnToadd], result, stringsAsFactors = F)

    ### list of permutated gene names
    genestopermut <- rownames(expMat)
    geneList <- vector("list", nPermu)
    names(geneList) <- paste0("permutation", 1:length(geneList))
    for(i in 1:length(geneList)){
      geneList[[i]] <- sample(genestopermut, replace = F)
    }

    #sample_pairInteraction <- sample(1:nrow(result), 1) ### scelgo una coppia a caso
    sample_pairInteraction <- 1 ### scelgo la prima coppia

    require(doMC)
    registerDoMC(ncore)
    ans <- foreach(gg = 1:length(geneList)) %dopar% {
      ddataA <- expMat[, cellCusterA]
      rownames(ddataA) <- geneList[[gg]]
      n <- round(nrow(ddataA) * topRank/100)
      ddataA_ <- apply(ddataA, 2, function(x) rank(-x))
      ddataA_ <- apply(ddataA_, 2, function(x) x <= n)
      ggenes <- result[sample_pairInteraction, c("partnerA1", "partnerA2", "partnerA3")]
      ggenes <- as.vector(ggenes[!is.na(ggenes)])
      tmp_genes <- ddataA_[ggenes, , drop = F]
      tmp <- colSums(tmp_genes)
      rnk_permuA <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)

      ddataB <- expMat[, cellCusterB]
      rownames(ddataB) <- geneList[[gg]]
      n <- round(nrow(ddataB) * topRank/100)
      ddataB_ <- apply(ddataB, 2, function(x) rank(-x))
      ddataB_ <- apply(ddataB_, 2, function(x) x <= n)
      ggenes <- result[sample_pairInteraction, c("partnerB1", "partnerB2", "partnerB3")]
      ggenes <- as.vector(ggenes[!is.na(ggenes)])
      tmp_genes <- ddataB_[ggenes, , drop = F]
      tmp <- colSums(tmp_genes)
      rnk_permuB <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)

      rnkInteraction <- mean(rnk_permuA, rnk_permuB)
      return(rnkInteraction)
    }

    permutatedPvalue <- unlist(ans)
    SCOREpValue <- c()
    for(k in 1:nrow(result)){
      SCOREpValue <- c(SCOREpValue, sum(permutatedPvalue > result[k, "SCORE"])/nPermu)
    }
    result <- data.frame(result, SCOREpValue, stringsAsFactors = F)
    result <- result[result$SCOREpValue <= pvalueCutoff, ]

    columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- result[, setdiff(colnames(result), columnToremove)]
    scTHIresult$result <- result
    save(scTHIresult, file = paste0(fileNameBase, "_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }

  message("Interaction pair:")
  print(rownames(scTHIresult$result))
  return(scTHIresult)
}


scTHI.plotResult <- function(scTHIresult, cexNames = 0.8, plotType = c("score", "pair"), legendPos = c("topleft", "topright")){

  #' scTHI.plotResults
  #'
  #' This function  ...
  #'
  #' @export
  #' @examples
  #'
  #' scTHI.plotResult

  result <- scTHIresult$result
  legendPos <- match.arg(legendPos)

  if(plotType == "score"){
    tmp <- as.matrix(result[, c("SCORE", "SCOREpValue")])
    pvalues <- tmp[,"SCOREpValue"]
    for(i in 1:length(pvalues)){
      if(tmp[i, "SCOREpValue"] > 0.05){
        pvalues[i] <- "ns"
      }
      if(tmp[i, "SCOREpValue"] <= 0.05){
        pvalues[i] <- "*"
      }
      if(tmp[i, "SCOREpValue"] <= 0.01){
        pvalues[i] <- "**"
      }
      if(tmp[i, "SCOREpValue"] <= 0.001){
        pvalues[i] <- "***"
      }
      if(tmp[i, "SCOREpValue"] <= 0.0001){
        pvalues[i] <- "****"
      }
    }

    barplot(t(tmp[, "SCORE"]), beside = T, ylim = c(0,1.1),  ylab = "scTHI Score", col = c("lightseagreen"), cex.names = cexNames)
    text(x = seq(1.5, by = 2, length.out = nrow(tmp)), y = 1.05, labels =  pvalues, cex = 2)

    #sizedot <- tmp[, "SCOREpValue"]
    #dotchart(tmp[, "SCORE"], xlim = c(0,1.1),  xlab = "scTHI Score", col = c("lightseagreen"), cex = 1, pch = 19, pt.cex = -log10(sizedot)*4, lcolor = "slategray3")
  }

  if(plotType == "pair"){
    tmp <- as.matrix(result[, grep("rnk", colnames(result))])
    barplot(t(tmp), beside = T, ylim = c(0,1.1),  ylab = "% of Cells", col = c("#F8766D", "#00BFC4"), cex.names = cexNames)
    legend(legendPos, legend = c("ClusterA", "ClusterB"), col = c("#F8766D", "#00BFC4"), cex = 1, pch = 15, box.lty=0)
  }
}


scTHI.runTsne <- function(scTHIresult){


  #' scTHI.runTsne
  #'
  #' This function  ...
  #'
  #' @export
  #' @examples
  #'
  #' scTHI.runTsne


  expMat <- scTHIresult$expMat

  ##create tsne
  eta <- .1
  filter <- apply(expMat, 1, function(x) sum(quantile(x, probs = c(1 - eta, eta)) * c(1, -1)))
  #fivenum(filter)
  #plot(density(log2(filter)))
  cutoff <- density(log2(filter))
  cutoff <- data.frame(x = cutoff$x, y = cutoff$y)
  cutoff <- cutoff$x[which.max(cutoff$y)]
  foldChange <- 2^cutoff
  #sum(filter > foldChange) # 4567
  variableGenes <- names(filter)[filter > foldChange]
  expMat <- expMat[variableGenes, ]

  require("Rtsne")
  expMatT <- t(expMat)
  set.seed(1) ### for reproducibility
  tsne_out <- Rtsne(expMatT)
  tsneData <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], Sample= colnames(expMat), stringsAsFactors = F)
  rownames(tsneData) <- colnames(expMat)
  scTHIresult <- c(list(tsneData), scTHIresult)
  names(scTHIresult)[1] <- "tsneData"
  return(scTHIresult)
}


scTHI.plotCluster <- function(scTHIresult, cexPoint = 0.8, legendPos = c("topleft", "topright", "bottomright", "bottomleft")){

  #' scTHI.plotCluster
  #'
  #' This function  ...
  #'
  #' @export
  #' @examples
  #'
  #' scTHI.plotCluster

  tsneData <- scTHIresult$tsneData
  legendPos <- match.arg(legendPos)

  nCluster <- length(scTHIresult) - 3
  ggplot_colors <- function(n){
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  Colors <- ggplot_colors(n = nCluster)
  names(Colors) <- names(scTHIresult)[4:length(scTHIresult)]
  ClusterColors <- rep("gray80", nrow(tsneData))
  names(ClusterColors) <- rownames(tsneData)
  for(i in 1:length(Colors)){
    ClusterColors[scTHIresult[[names(Colors)[i]]]] <- Colors[i]
  }
  plot(tsneData[,1:2], pch = 16,  cex = cexPoint , col = ClusterColors[rownames(tsneData)], main = "", xlab = "tSNE 1", ylab = "tSNE 2")
  legend(legendPos, inset=.02, legend = names(Colors), col = Colors, cex = 1, pch = 15, box.lty=0)
}


scTHI.plotPairs <- function(scTHIresult, cexPoint = 0.8, interactionToplot){

  #' scTHI.plotPairs
  #'
  #' This function  ...
  #'
  #' @export
  #' @examples
  #'
  #' scTHI.plotPairs


  tsneData <- scTHIresult$tsneData
  result <- scTHIresult$result
  expMat <- scTHIresult$expMat

  genesToplotA <- result[interactionToplot, c("partnerA")]
  genesToplotA <- unlist(strsplit(genesToplotA, ":"))
  genesToplotB <- result[interactionToplot, c("partnerB")]
  genesToplotB <- unlist(strsplit(genesToplotB, ":"))



  #################  create list colori
  colByValue_ <- function (x, col, range = NA, breaks = NA, cex.axis = 2, las = 1,
                           ...)
  {
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
      ff <- seq(min(x, na.rm = T), max(x, na.rm = T), length = length(col))
      bg2 <- apply(as.matrix(as.numeric(unlist(x))), 1, function(x) rank(c(ff,
                                                                           x), ties.method = "min")[length(col) + 1])
      dens <- matrix(bg2, nrow(x), ncol(x))
      result <- matrix(col[dens], nrow = nrow(x), ncol = ncol(x))
      row.names(result) <- row.names(x)
      # image(x = 1:2, y = as.matrix(ff), z = t(ff), col = col,
      #       xaxt = "n", ylab = "", las = las, xlab = "", xlim = c(1,
      #                                                             4), bty = "n", ...)
      return(result)
    }
    else {
      temp <- cut(as.numeric(unlist(x)), breaks = breaks, include.lowest = T)
      if (length(col) != length(levels(temp))) {
        stop("length:col != length: cut result")
      }
      result <- matrix(col[as.numeric(temp)], nrow = nrow(x),
                       ncol = ncol(x))
      row.names(result) <- row.names(x)
      # image(x = 1:2, y = as.matrix(1:(length(breaks) - 1)),
      #       z = t(1:(length(breaks) - 1)), col = col, xaxt = "n",
      #       yaxt = "n", ylab = "", xlab = "", xlim = c(0, 3),
      #       ...)
      # axis(2, at = 1:(length(breaks) - 1), labels = levels(temp),
      #      las = las, cex.axis = cex.axis)
      return(result)
    }
  }

  list_colori_A <- vector("list", length(genesToplotA))
  names(list_colori_A) <- genesToplotA
  for(i in 1:length(genesToplotA)){
    tmp_exp <- t(expMat[genesToplotA[i], colnames(expMat), drop = F])
    mean_exp_samples <- mean(tmp_exp)

    ##### range gene x gene
    range_exp <- range(tmp_exp)
    bre <- seq(min(tmp_exp[, genesToplotA[i]]) - 0.01, max(tmp_exp[, genesToplotA[i]]) + 0.01, 0.01)
    colori <- colByValue_(tmp_exp, col= colorRampPalette(c('gray80', 'coral','red'))(length(bre)-1), breaks= bre, cex.axis=0.8)
    colori <- colori[rownames(tsneData), , drop = F]
    tmp <- list(colori, range_exp, mean_exp_samples)
    names(tmp) <- c("colori", "range_exp", "mean_exp_samples")
    list_colori_A[[i]] <- tmp
  }

  list_colori_B <- vector("list", length(genesToplotB))
  names(list_colori_B) <- genesToplotB
  for(i in 1:length(genesToplotB)){
    tmp_exp <- t(expMat[genesToplotB[i], colnames(expMat), drop = F])
    mean_exp_samples <- mean(tmp_exp)

    ##### range gene x gene
    range_exp <- range(tmp_exp)
    bre <- seq(min(tmp_exp[, genesToplotB[i]]) - 0.01, max(tmp_exp[, genesToplotB[i]]) + 0.01, 0.01)
    colori <- colByValue_(tmp_exp, col= colorRampPalette(c('gray80', 'coral','red'))(length(bre)-1), breaks= bre, cex.axis=0.8)
    colori <- colori[rownames(tsneData), , drop = F]
    tmp <- list(colori, range_exp, mean_exp_samples)
    names(tmp) <- c("colori", "range_exp", "mean_exp_samples")
    list_colori_B[[i]] <- tmp
  }

  ################### plot
  par(mfrow= c(1, length(c(genesToplotA, genesToplotB))), mar=c(6, 4, 3, 2), cex.main = 1.5, cex.sub = 1.2, col.main = "black", col.sub = "gray30", font.main = 2, font.sub = 3)
  for(i in 1:length(list_colori_A)){
    colori <- list_colori_A[[i]]$colori
    range_exp <- round(list_colori_A[[i]]$range_exp, digits = 2)
    mean_exp_samples <- round(list_colori_A[[i]]$mean_exp_samples, digits = 2)
    plot(tsneData[,1:2], pch = 16,  cex = cexPoint , col = colori,
         main = genesToplotA[i], xlab = "", ylab = "",
         sub = paste0("Range: ", range_exp[1], " to ",  range_exp[2],", MeanExp: ", mean_exp_samples))
  }

  for(i in 1:length(list_colori_B)){
    colori <- list_colori_B[[i]]$colori
    range_exp <- round(list_colori_B[[i]]$range_exp, digits = 2)
    mean_exp_samples <- round(list_colori_B[[i]]$mean_exp_samples, digits = 2)
    plot(tsneData[,1:2], pch = 16,  cex = cexPoint , col = colori,
         main = genesToplotB[i], xlab = "", ylab = "",
         sub = paste0("Range: ", range_exp[1], " to ",  range_exp[2],", MeanExp: ", mean_exp_samples))

  }
}


TME.classification <- function(expMat, ncore = 48, minLenGeneSet = 10, pvalFilter = FALSE, fdrFilter = TRUE, pvalCutoff = 0.01, nesCutoff = 0.58, nNES = 1){

  #' TME.classification
  #'
  #' This function  ...
  #'
  #' @export
  #' @examples
  #'
  #' TME.classification


  library(yaGST)
  means <- rowMeans(expMat)
  sds <- apply(expMat, 1, sd)
  library(doMC)
  registerDoMC(ncore)
  message("Computing ssMWW-GST ...")
  ans <- foreach(ss = 1:ncol(expMat)) %dopar% {
    currentSample <- (expMat[, ss] - means)/sds
    rankedList <- sort(currentSample, decreasing = T)
    aMwwGST <- lapply(signatures, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = minLenGeneSet, alternative = "two.sided"))
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
    ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
    return(ans)
  }

  NES <- sapply(ans, function(x) x$tmp_NES)
  pValue <- sapply(ans, function(x) x$tmp_pValue)
  colnames(NES) <- colnames(pValue) <- colnames(expMat)
  FDR <- apply(pValue, 2, function(x) p.adjust(x, method = "fdr"))

  if(pvalFilter == TRUE){
    NES[pValue >= pvalCutoff] <- 0
  }
  if(fdrFilter == TRUE){
    NES[FDR >= pvalCutoff] <- 0
  }
  NES[NES < nesCutoff] <- 0
  if(nNES == 1){
    Class <- apply(NES, 2, function(x) names(which.max(x)))
  }
  if(nNES != 1){
    Class <- apply(NES, 2, function(x) names(sort(x, decreasing = T))[nNES])
  }
  Class[colSums(NES != 0) < nNES] <- "nc"

  phenotype <- signaturesColors[Class, ]
  rownames(phenotype) <- names(Class)

  Class <- phenotype$ALLPhenotypeFinal
  names(Class) <- rownames(phenotype)

  ClassLegend <- phenotype$Color
  names(ClassLegend) <- phenotype$ALLPhenotypeFinal
  ClassLegend <- ClassLegend[!duplicated(ClassLegend)]

  print(sort(table(Class), decreasing = T))
  Classification <- list(Class,  ClassLegend)
  names(Classification) <- c("Class", "ClassLegend")
  return(Classification)
}


TME.plot <- function(tsneData, Class, cexPoint = 0.8){
  #' TME.plot
  #'
  #' This function  ...
  #'
  #' @export
  #' @examples
  #'
  #' TME.plot

  ClassColor <- Class$ClassLegend[Class$Class]
  names(ClassColor) <- names(Class$Class)

  TmeColors <- rep("gray80", nrow(tsneData))
  names(TmeColors) <- rownames(tsneData)
  TmeColors[names(ClassColor)] <- ClassColor
  par(oma=c(0, 0, 0, 6))
  plot(tsneData[,1:2], pch = 16,  cex = cexPoint , col = TmeColors[rownames(tsneData)], main = "", xlab = "tSNE 1", ylab = "tSNE 2")
  legend(par('usr')[2], par('usr')[4],  bty='n', xpd=NA, legend = names(Class$ClassLegend), col = Class$ClassLegend, cex = 0.9, pch = 15, box.lty=0)
}
