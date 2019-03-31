
scTHI <- function(expMat, cellCusterA, cellCusterB, cellCusterAName, cellCusterBName, topRank = 10, fileNameBase = "CustomCellPhone", filterCutoff, 
                                 foldChange = TRUE, useData = FALSE, dataToUse, cellsToCompare, cellsToCompareNAME, expMatToCompare,
                                 PValue = TRUE, nPermu = 1000, ncore = 8){
  
  #' scTHI
  #'
  #' This function allows you to express your love of cats.
  #' @param expMat expression matrix.
  #' @keywords interaction
  #' @export
  #' @examples
  #' scTHI
  
  ####################### check rownames expMat
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
      print(i)
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
    save(result, file = paste0(fileNameBase, "result_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }
  
  ######################### compute permutations ##########
  if(PValue == TRUE){
    message("Computing permutation....")
    
    if(useData == TRUE){
      result <- dataToUse
    }
    if(nrow(result) == 0) stop("No result!")
    
    load("~/result/data_interaction/interaction_tableComplete.RData")
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
    result <- result[result$SCOREpValue < 0.05, ]
    
    columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- result[, setdiff(colnames(result), columnToremove)]
    save(result, file = paste0(fileNameBase, "result_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }
  
  ######################### compute folde change ##########
  if(foldChange == TRUE){
    message("Computing fold change....")
    
    if(useData == TRUE){
      result <- dataToUse
    }
    if(nrow(result) == 0) stop("No result!")
    
    load("~/result/data_interaction/interaction_tableComplete.RData")
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
        log2FC_partnerB[i] <- mean(rowMeans(expMat[ggenesB, cellCusterA, drop = F])) - mean(rowMeans(expMatToCompare[ggenesB, cellsToCompare, drop = F]))
      }
      
    }
    FC <- cbind(log2FC_partnerA, log2FC_partnerB)
    colnames(FC) <- paste0(colnames(FC), "__", cellsToCompareNAME)
    result <- data.frame(result, FC, stringsAsFactors = F)
    
    columnToremove <- c("partnerA1", "partnerA2", "partnerA3", "partnerB1", "partnerB2", "partnerB3")
    result <- result[, setdiff(colnames(result), columnToremove)]
    save(result, file = paste0(fileNameBase, "result_", cellCusterAName, "VS", cellCusterBName, ".RData"))
  }
  
  ###############
  return(result)
}