
#' scTHI <- function(expMat, cellCusterA, cellCusterB, cellCusterAName, cellCusterBName, topRank = 10, fileNameBase = "CustomCellPhone", filterCutoff,
#'                                  foldChange = TRUE, useData = FALSE, dataToUse, cellsToCompare, cellsToCompareNAME, expMatToCompare,
#'                                  PValue = TRUE, nPermu = 1000, ncore = 8){
#'
#'   #' scTHI
#'   #'
#'   #' This function computes a set of ligand-receptor interactions from a single cell gene expression matrix.
#'   #' you must specify at least two clusters of cells (for example tumor cells and immune cells)
#'   #' @param expMat expression matrix.
#'   #' @param clusterA vector of columns of expMat that belong to the first cluster
#'   #' @param clusterB vector of columns of expMat that belong to the second cluster
#'   #' @param cellCusterBName
#'   #' @keywords interaction
#'   #' @export
#'   #' @examples
#'   #'
#'   #' scTHI
#'
#'   ####################### check rownames expMat
#'   tmp_check <- interaction_table[, c("partnerA1" ,"partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")]
#'   tmp_check <- unique(tmp_check[!is.na(tmp_check)])
#'   if(sum(tmp_check %in% rownames(expMat)) != length(tmp_check)){
#'     tmp_genes <- setdiff(tmp_check, rownames(expMat))
#'     tmp <- which(interaction_table$partnerA1 %in% tmp_genes | interaction_table$partnerA2 %in% tmp_genes | interaction_table$partnerA3 %in% tmp_genes | interaction_table$partnerB1 %in% tmp_genes | interaction_table$partnerB2 %in% tmp_genes | interaction_table$partnerB3 %in% tmp_genes)
#'     interaction_table <- interaction_table[-tmp, ]
#'     message("not all interaction genes are in expMat")
#'   }
#'   if(nrow(interaction_table) == 0) stop("No interaction genes to test")
#'
#'   ######################### compute score ##################
#'   if(useData == FALSE){
#'     message("Computing score....")
#'
#'     ddataA <- expMat[, cellCusterA]
#'     n <- round(nrow(ddataA) * topRank/100)
#'     ddataA_ <- apply(ddataA, 2, function(x) rank(-x))
#'     ddataA_ <- apply(ddataA_, 2, function(x) x <= n)
#'     rnkA <- rep(0, nrow(interaction_table))
#'     expValueA <- rep(0, nrow(interaction_table))
#'     names(rnkA) <- names(expValueA) <- rownames(interaction_table)
#'     for(i in 1:nrow(interaction_table)){
#'       print(i)
#'       ggenes <- interaction_table[i, c("partnerA1", "partnerA2", "partnerA3")]
#'       ggenes <- as.vector(ggenes[!is.na(ggenes)])
#'       tmp_genes <- ddataA_[ggenes, , drop = F]
#'       tmp <- colSums(tmp_genes)
#'       rnkA[i] <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
#'       expValueA[i] <- paste(round(rowMeans(ddataA[ggenes, cellCusterA, drop = F]), digits = 2), collapse = " # ")
#'     }
#'     resultA <- data.frame(rnkPartnerA = rnkA, expValueA = expValueA, stringsAsFactors = F)
#'     colnames(resultA) <- paste0(colnames(resultA), "_", cellCusterAName)
#'
#'     ddataB <- expMat[, cellCusterB]
#'     n <- round(nrow(ddataB) * topRank/100)
#'     ddataB_ <- apply(ddataB, 2, function(x) rank(-x))
#'     ddataB_ <- apply(ddataB_, 2, function(x) x <= n)
#'     rnkB <- rep(0, nrow(interaction_table))
#'     expValueB <- rep(0, nrow(interaction_table))
#'     names(rnkB) <- names(expValueB) <- rownames(interaction_table)
#'     for(i in 1:nrow(interaction_table)){
#'       print(i)
#'       ggenes <- interaction_table[i, c("partnerB1", "partnerB2", "partnerB3")]
#'       ggenes <- as.vector(ggenes[!is.na(ggenes)])
#'       tmp_genes <- ddataB_[ggenes, , drop = F]
#'       tmp <- colSums(tmp_genes)
#'       rnkB[i] <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
#'       expValueB[i] <- paste(round(rowMeans(ddataB[ggenes, cellCusterB, drop = F]), digits = 2), collapse = " # ")
#'     }
#'     resultB <- data.frame(rnkPartnerB = rnkB, expValueB = expValueB, stringsAsFactors = F)
#'     colnames(resultB) <- paste0(colnames(resultB), "_", cellCusterBName)
#'
#'     SCORE <- rowMeans(cbind(rnkA, rnkB))
#'     result <- data.frame(interaction_table, resultA, resultB, SCORE)
#'
#'     #### remove low rank-score interaction
#'     interestColumn <- c(paste0("rnkPartnerA_", cellCusterAName), paste0("rnkPartnerB_", cellCusterBName))
#'     tmp <- rowSums(result[, interestColumn] < filterCutoff)
#'     result[tmp != 0, "SCORE"] <- NA
#'     result <- result[!is.na(result$SCORE), ]
#'
#'     #### prepare final result
#'     result <- result[order(result$SCORE, decreasing = T), ]
#'     result <- data.frame(interationPair = rownames(result), result, stringsAsFactors = F)
#'     columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
#'     result <- result[, setdiff(colnames(result), columnToremove)]
#'
#'     if(nrow(result) == 0) stop("No result!")
#'     save(result, file = paste0(fileNameBase, "result_", cellCusterAName, "VS", cellCusterBName, ".RData"))
#'   }
#'
#'   ######################### compute permutations ##########
#'   if(PValue == TRUE){
#'     message("Computing permutation....")
#'
#'     if(useData == TRUE){
#'       result <- dataToUse
#'     }
#'     if(nrow(result) == 0) stop("No result!")
#'
#'     #load("~/result/data_interaction/interaction_tableComplete.RData")
#'     interaction_table <- interaction_table[rownames(result), ]
#'     columnToadd <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
#'     result <- data.frame(interaction_table[, columnToadd], result, stringsAsFactors = F)
#'
#'     ### list of permutated gene names
#'     genestopermut <- rownames(expMat)
#'     geneList <- vector("list", nPermu)
#'     names(geneList) <- paste0("permutation", 1:length(geneList))
#'     for(i in 1:length(geneList)){
#'       geneList[[i]] <- sample(genestopermut, replace = F)
#'     }
#'
#'     #sample_pairInteraction <- sample(1:nrow(result), 1) ### scelgo una coppia a caso
#'     sample_pairInteraction <- 1 ### scelgo la prima coppia
#'
#'     require(doMC)
#'     registerDoMC(ncore)
#'     ans <- foreach(gg = 1:length(geneList)) %dopar% {
#'       ddataA <- expMat[, cellCusterA]
#'       rownames(ddataA) <- geneList[[gg]]
#'       n <- round(nrow(ddataA) * topRank/100)
#'       ddataA_ <- apply(ddataA, 2, function(x) rank(-x))
#'       ddataA_ <- apply(ddataA_, 2, function(x) x <= n)
#'       ggenes <- result[sample_pairInteraction, c("partnerA1", "partnerA2", "partnerA3")]
#'       ggenes <- as.vector(ggenes[!is.na(ggenes)])
#'       tmp_genes <- ddataA_[ggenes, , drop = F]
#'       tmp <- colSums(tmp_genes)
#'       rnk_permuA <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
#'
#'       ddataB <- expMat[, cellCusterB]
#'       rownames(ddataB) <- geneList[[gg]]
#'       n <- round(nrow(ddataB) * topRank/100)
#'       ddataB_ <- apply(ddataB, 2, function(x) rank(-x))
#'       ddataB_ <- apply(ddataB_, 2, function(x) x <= n)
#'       ggenes <- result[sample_pairInteraction, c("partnerB1", "partnerB2", "partnerB3")]
#'       ggenes <- as.vector(ggenes[!is.na(ggenes)])
#'       tmp_genes <- ddataB_[ggenes, , drop = F]
#'       tmp <- colSums(tmp_genes)
#'       rnk_permuB <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
#'
#'       rnkInteraction <- mean(rnk_permuA, rnk_permuB)
#'       return(rnkInteraction)
#'     }
#'
#'     permutatedPvalue <- unlist(ans)
#'     SCOREpValue <- c()
#'     for(k in 1:nrow(result)){
#'       SCOREpValue <- c(SCOREpValue, sum(permutatedPvalue > result[k, "SCORE"])/nPermu)
#'     }
#'     result <- data.frame(result, SCOREpValue, stringsAsFactors = F)
#'     result <- result[result$SCOREpValue < 0.05, ]
#'
#'     columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
#'     result <- result[, setdiff(colnames(result), columnToremove)]
#'     save(result, file = paste0(fileNameBase, "result_", cellCusterAName, "VS", cellCusterBName, ".RData"))
#'   }
#'
#'   ######################### compute folde change ##########
#'   if(foldChange == TRUE){
#'     message("Computing fold change....")
#'
#'     if(useData == TRUE){
#'       result <- dataToUse
#'     }
#'     if(nrow(result) == 0) stop("No result!")
#'
#'     #load("~/result/data_interaction/interaction_tableComplete.RData")
#'     interaction_table <- interaction_table[rownames(result), ]
#'     columnToadd <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
#'     result <- data.frame(interaction_table[, columnToadd], result, stringsAsFactors = F)
#'
#'     ########### check rownames expMatToCompare
#'     tmp_check <- result[, c("partnerA1" ,"partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")]
#'     tmp_check <- unique(tmp_check[!is.na(tmp_check)])
#'     if(sum(tmp_check %in% rownames(expMatToCompare)) != length(tmp_check)){
#'       message("Not all interaction genes are in expMatToCompare")
#'     }
#'
#'     log2FC_partnerA <- rep(NA, nrow(result))
#'     log2FC_partnerB <- rep(NA, nrow(result))
#'     for(i in 1:nrow(result)){
#'       ggenesA <- result[i, c("partnerA1", "partnerA2" , "partnerA3")]
#'       ggenesA <- as.vector(ggenesA[!is.na(ggenesA)])
#'
#'       ggenesB <- result[i, c("partnerB1", "partnerB2" , "partnerB3")]
#'       ggenesB <- as.vector(ggenesB[!is.na(ggenesB)])
#'
#'       if(sum(ggenesA %in% rownames(expMatToCompare)) == length(ggenesA)){
#'         log2FC_partnerA[i] <- mean(rowMeans(expMat[ggenesA, cellCusterA, drop = F])) - mean(rowMeans(expMatToCompare[ggenesA, cellsToCompare, drop = F]))
#'       }
#'       if(sum(ggenesB %in% rownames(expMatToCompare)) == length(ggenesB)){
#'         log2FC_partnerB[i] <- mean(rowMeans(expMat[ggenesB, cellCusterA, drop = F])) - mean(rowMeans(expMatToCompare[ggenesB, cellsToCompare, drop = F]))
#'       }
#'
#'     }
#'     FC <- cbind(log2FC_partnerA, log2FC_partnerB)
#'     colnames(FC) <- paste0(colnames(FC), "__", cellsToCompareNAME)
#'     result <- data.frame(result, FC, stringsAsFactors = F)
#'
#'     columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
#'     result <- result[, setdiff(colnames(result), columnToremove)]
#'     save(result, file = paste0(fileNameBase, "result_", cellCusterAName, "VS", cellCusterBName, ".RData"))
#'   }
#'
#'   ###############
#'   return(result)
#' }
#'
#' scTHI <- function(expMat, cellCusterA, cellCusterB, cellCusterAName, cellCusterBName, topRank = 10, fileNameBase = "CustomCellPhone", filterCutoff,
foldChange = TRUE, useData = FALSE, dataToUse, cellsToCompare, cellsToCompareNAME, expMatToCompare,
PValue = TRUE, nPermu = 1000, ncore = 8){

  #' scTHI
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
  #load("/storage/gluster/vol1/SHARED/HOMEFOLDERS/caruso/scProject/NYnontumor/Phone/PhoneCustom/package_data/interaction_tableComplete_1144.RData") #### da rimuovere x il pacchetto
  tmp_check <- interaction_table[, c("partnerA1" ,"partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")]
  tmp_check <- unique(tmp_check[!is.na(tmp_check)])
  if(sum(tmp_check %in% rownames(expMat)) != length(tmp_check)){
    tmp_genes <- setdiff(tmp_check, rownames(expMat))
    tmp <- which(interaction_table$partnerA1 %in% tmp_genes | interaction_table$partnerA2 %in% tmp_genes | interaction_table$partnerA3 %in% tmp_genes | interaction_table$partnerB1 %in% tmp_genes | interaction_table$partnerB2 %in% tmp_genes | interaction_table$partnerB3 %in% tmp_genes)
    interaction_table <- interaction_table[-tmp, ]
    message("not all interaction genes are in expMat")
  }
  if(nrow(interaction_table) == 0) stop("No interaction genes to test")

  ######################### compute score ##################
  if(useData == FALSE){
    message("Computing score....")

    ddataA <- expMat[, cellCusterA]
    n <- round(nrow(ddataA) * topRank/100)
    ddataA_ <- apply(ddataA, 2, function(x) rank(-x))
    ddataA_ <- apply(ddataA_, 2, function(x) x <= n)
    rnkA <- rep(0, nrow(interaction_table))
    expValueA <- rep(0, nrow(interaction_table))
    names(rnkA) <- names(expValueA) <- rownames(interaction_table)
    for(i in 1:nrow(interaction_table)){
      #print(i)
      ggenes <- interaction_table[i, c("partnerA1", "partnerA2", "partnerA3")]
      ggenes <- as.vector(ggenes[!is.na(ggenes)])
      tmp_genes <- ddataA_[ggenes, , drop = F]
      tmp <- colSums(tmp_genes)
      rnkA[i] <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
      expValueA[i] <- paste(round(rowMeans(ddataA[ggenes, cellCusterA, drop = F]), digits = 2), collapse = " # ")
    }
    resultA <- data.frame(rnkPartnerA = rnkA, expValueA = expValueA, stringsAsFactors = F)
    colnames(resultA) <- paste0(colnames(resultA), "_", cellCusterAName)

    ddataB <- expMat[, cellCusterB]
    n <- round(nrow(ddataB) * topRank/100)
    ddataB_ <- apply(ddataB, 2, function(x) rank(-x))
    ddataB_ <- apply(ddataB_, 2, function(x) x <= n)
    rnkB <- rep(0, nrow(interaction_table))
    expValueB <- rep(0, nrow(interaction_table))
    names(rnkB) <- names(expValueB) <- rownames(interaction_table)
    for(i in 1:nrow(interaction_table)){
      print(i)
      ggenes <- interaction_table[i, c("partnerB1", "partnerB2", "partnerB3")]
      ggenes <- as.vector(ggenes[!is.na(ggenes)])
      tmp_genes <- ddataB_[ggenes, , drop = F]
      tmp <- colSums(tmp_genes)
      rnkB[i] <- round(sum(tmp == nrow(tmp_genes))/ncol(tmp_genes), digits = 2)
      expValueB[i] <- paste(round(rowMeans(ddataB[ggenes, cellCusterB, drop = F]), digits = 2), collapse = " # ")
    }
    resultB <- data.frame(rnkPartnerB = rnkB, expValueB = expValueB, stringsAsFactors = F)
    colnames(resultB) <- paste0(colnames(resultB), "_", cellCusterBName)

    SCORE <- rowMeans(cbind(rnkA, rnkB))
    result <- data.frame(interaction_table, resultA, resultB, SCORE)

    #### remove low rank-score interaction
    interestColumn <- c(paste0("rnkPartnerA_", cellCusterAName), paste0("rnkPartnerB_", cellCusterBName))
    tmp <- rowSums(result[, interestColumn] < filterCutoff)
    result[tmp != 0, "SCORE"] <- NA
    result <- result[!is.na(result$SCORE), ]

    #### prepare final result
    result <- result[order(result$SCORE, decreasing = T), ]
    result <- data.frame(interationPair = rownames(result), result, stringsAsFactors = F)
    columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- result[, setdiff(colnames(result), columnToremove)]

    if(nrow(result) == 0) stop("No result!")
    save(result, file = paste0(fileNameBase, "_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }

  ######################### compute permutations ##########
  if(PValue == TRUE){
    message("Computing permutation....")

    if(useData == TRUE){
      result <- dataToUse
    }
    if(nrow(result) == 0) stop("No result!")

    #load("/storage/gluster/vol1/SHARED/HOMEFOLDERS/caruso/scProject/NYnontumor/Phone/PhoneCustom/package_data/interaction_tableComplete_1144.RData") #### da rimuovere per il pacchetto
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

    ####sample_pairInteraction <- sample(1:nrow(result), 1) ### scelgo una coppia a caso
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
    result <- result[result$SCOREpValue < 0.05, ]

    columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- result[, setdiff(colnames(result), columnToremove)]
    save(result, file = paste0(fileNameBase, "_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }

  ######################### compute folde change ##########
  if(foldChange == TRUE){
    message("Computing fold change....")

    if(useData == TRUE){
      result <- dataToUse
    }
    if(nrow(result) == 0) stop("No result!")

    #load("/storage/gluster/vol1/SHARED/HOMEFOLDERS/caruso/scProject/NYnontumor/Phone/PhoneCustom/package_data/interaction_tableComplete_1144.RData") #### da rimuovere per il pacchetto
    interaction_table <- interaction_table[rownames(result), ]
    columnToadd <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- data.frame(interaction_table[, columnToadd], result, stringsAsFactors = F)

    ########### check rownames expMatToCompare
    tmp_check <- result[, c("partnerA1" ,"partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")]
    tmp_check <- unique(tmp_check[!is.na(tmp_check)])
    if(sum(tmp_check %in% rownames(expMatToCompare)) != length(tmp_check)){
      message("Not all interaction genes are in expMatToCompare")
    }

    log2FC_partnerA <- rep(NA, nrow(result))
    log2FC_partnerB <- rep(NA, nrow(result))
    for(i in 1:nrow(result)){
      ggenesA <- result[i, c("partnerA1", "partnerA2" , "partnerA3")]
      ggenesA <- as.vector(ggenesA[!is.na(ggenesA)])

      ggenesB <- result[i, c("partnerB1", "partnerB2" , "partnerB3")]
      ggenesB <- as.vector(ggenesB[!is.na(ggenesB)])

      if(sum(ggenesA %in% rownames(expMatToCompare)) == length(ggenesA)){
        log2FC_partnerA[i] <- mean(rowMeans(expMat[ggenesA, cellCusterA, drop = F])) - mean(rowMeans(expMatToCompare[ggenesA, cellsToCompare, drop = F]))
      }
      if(sum(ggenesB %in% rownames(expMatToCompare)) == length(ggenesB)){
        log2FC_partnerB[i] <- mean(rowMeans(expMat[ggenesB, cellCusterB, drop = F])) - mean(rowMeans(expMatToCompare[ggenesB, cellsToCompare, drop = F]))
      }

    }
    FC <- cbind(log2FC_partnerA, log2FC_partnerB)
    colnames(FC) <- paste0(colnames(FC), "__", c(cellCusterAName, cellCusterBName), "VS",cellsToCompareNAME)
    result <- data.frame(result, FC, stringsAsFactors = F)

    columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- result[, setdiff(colnames(result), columnToremove)]
    save(result, file = paste0(fileNameBase, "_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }

  ###############
  return(result)
}
