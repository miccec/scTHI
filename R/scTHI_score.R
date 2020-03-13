checkInteractions <- function(expMat){
  tmp_check <-
    interaction_table[, c(
      "partnerA1",
      "partnerA2",
      "partnerA3",
      "partnerB1",
      "partnerB2",
      "partnerB3"
    )]
  tmp_check <- unique(tmp_check[!is.na(tmp_check)])
  if (sum(tmp_check %in% rownames(expMat)) != length(tmp_check)) {
    tmp_genes <- setdiff(tmp_check, rownames(expMat))
    tmp <- which(
      interaction_table$partnerA1 %in% tmp_genes |
        interaction_table$partnerA2 %in% tmp_genes |
        interaction_table$partnerA3 %in% tmp_genes |
        interaction_table$partnerB1 %in% tmp_genes |
        interaction_table$partnerB2 %in% tmp_genes |
        interaction_table$partnerB3 %in% tmp_genes
    )
    interaction_table <- interaction_table[-tmp, ]
    message("Warning: Not all interaction genes are in expMat")
  }
  if (nrow(interaction_table) == 0) {
    stop("ERROR: No interaction genes to test")
  }
  return(interaction_table)
}

rankingValues <- function(x, topRank = 10){
  n <- round(nrow(x) * topRank / 100)
  ddata_ <- apply(x, 2, function(x) {
    rank(-x)
  })
  ddata_ <- apply(ddata_, 2, function(x) {
    x <= n
  })
  return(ddata_)
}

getScore <- function(expMat, interaction_table, cellCuster, 
                     genes1columns, genes2columns,autocrineEffect){
  ddata <- expMat[, cellCuster]
  ddata_ <- rankingValues(ddata)
  
  rnk <- rep(0, nrow(interaction_table))
  expValue <- rep(0, nrow(interaction_table))
  names(rnk) <- rownames(interaction_table)
  names(expValue) <- rownames(interaction_table)
  
  ans <-  apply(interaction_table, 1, function(x){
    ggenes1 <- x[genes1columns]
    ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
    tmp_genes1 <- ddata_[ggenes1, , drop = FALSE]
    tmp1 <- colSums(tmp_genes1 == TRUE)
    
    if (autocrineEffect == FALSE) {
      ggenes2 <-
        x[genes2columns]
      ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
      tmp_genes2 <- ddata_[ggenes2, , drop = FALSE]
      tmp2 <- colSums(tmp_genes2 == FALSE)
      tmp <- tmp1 + tmp2
      
      rnk<- round(sum(tmp == (nrow(tmp_genes1) +
                                nrow(tmp_genes2))) / ncol(tmp_genes1), digits = 2)
      
      expValue <-
        paste(round(rowMeans(ddata[ggenes1, cellCuster,
                                   drop = FALSE
                                   ]), digits = 2), collapse = " # ")
      
    }
    if (autocrineEffect == TRUE)  {
      tmp <- tmp1
      rnk <-
        round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
              digits = 2
        )
      expValue <-
        paste(round(rowMeans(ddata[ggenes1, cellCuster,
                                   drop = FALSE
                                   ]), digits = 2), collapse = " # ")
    }
    list(rnk,expValue)
  })
  rresult <- data.frame(
    rnk = unlist(lapply(ans, function(x) x[[1]])),
    expValue = unlist(lapply(ans, function(x) x[[2]]))
  )
  
  return(rresult)
}


#' scTHI_score
#'
#' This function allows the user to compute a score for a set of
#' ligand-receptor pairs, from a single cell gene expression matrix,
#' and detect specific Tumor-Host interactions. You must specify at
#' least two clusters of cells (for example tumor cells and immune
#' cells).
#' @param expMat ScRNA-seq gene expression matrix where rows are genes
#'   presented with Hugo Symbols and columns are cells. Gene expression
#'   values should be counts or normalized counts.
#' @param cellCusterA Vector of columns of expMat that belong to the
#' first cluster.
#' @param cellCusterB Vector of columns of expMat that belong to the
#' second cluster.
#' @param cellCusterAName A character string labeling the clusterA.
#' @param cellCusterBName A character string labeling the clusterB.
#' @param topRank Filter threshold. Set to 10 (default) means that
#' each gene of the interaction pair will be considered as expressed
#' in a cell if it's in the top rank 10 percent.
#' @param autocrineEffect if TRUE remove the paracrine filter
#' @param fileNameBase Project name.
#' @param filterCutoff Score threshold (default is 0.50). For each
#'   interaction pair, if the score calculated (for the partnerA
#'   or partnerB)
#'   will be less than filterCutoff the interaction pair will be
#'   discarded.
#' @param PValue Logical, set to TRUE (default) compute statistical
#'   iterations. If p.value < 0.05, the value will be returned.
#' @param pvalueCutoff cutoff of the p-value
#' @param nPermu Number of iterations to perform (default is 1000).
#' @param ncore Number of processors to use.
#' @keywords interaction
#' @examples
#'
#' ####################### example of scTHI_score
#' library(scTHI.data)
#' data(scExample)
#' result <-  scTHI_score(scExample,
#'       cellCusterA = colnames(scExample)[1:30],
#'       cellCusterB = colnames(scExample)[31:100],
#'       cellCusterAName = "ClusterA",
#'       cellCusterBName = "ClusterB", filterCutoff = 0,
#'      pvalueCutoff = 1, nPermu = 100, ncore = 8)
#'
#' @return A list of results, with four items: result (data.frame),
#'   expMat (matrix), clusterA (character),  clusterA (character)
#' @export
#' scTHI_score

scTHI_score <-
  function(expMat,
           cellCusterA,
           cellCusterB,
           cellCusterAName,
           cellCusterBName,
           topRank = 10,
           autocrineEffect = TRUE,
           fileNameBase = "scTHI",
           filterCutoff = 0.5,
           PValue = TRUE,
           pvalueCutoff = 0.05,
           nPermu = 1000,
           ncore = 8) {
    
   
    ############ check all interaction genes in expMat ############
    interaction_table <- checkInteractions(expMat)
    
    #################### scTHI score ##############################
    message(paste(
      "Computing score for",
      nrow(interaction_table),
      "interaction pairs..."
    ))
    
    ##### ddataA
    resultA <- getScore(expMat, interaction_table = interaction_table,
                        cellCuster = cellCusterA,
                        genes1columns= c("partnerA1", "partnerA2", 
                                         "partnerA3"),
                        genes2columns=c("partnerB1", "partnerB2", 
                                        "partnerB3"),
                        autocrineEffect)
    colnames(resultA) <-
      paste0(c("rnkPartnerA", "expValueA"), "_", cellCusterAName)
    print(paste("Computed ranked values for partner A"))
    
    ###### ddataB
    resultB <- getScore(expMat, interaction_table = interaction_table,
                        cellCuster = cellCusterB,
                        genes1columns= c("partnerB1", "partnerB2", 
                                         "partnerB3"),
                        genes2columns=c("partnerA1", "partnerA2", 
                                        "partnerA3"),
                        autocrineEffect)
    colnames(resultB) <-
      paste0(c("rnkPartnerB", "expValueB"), "_", cellCusterBName)
    print(paste("Computed ranked values for partner B"))
    
    
    #################### merge results #################
    SCORE <- rowMeans(cbind(resultA$rnkPartnerA, resultB$rnkPartnerB))
    result <- data.frame(interaction_table, resultA, resultB, SCORE)
    
    ############### filter low-score int ###############
    message("Removing low score interactions...")
    interestColumn <- c(
      paste0("rnkPartnerA_", cellCusterAName),
      paste0("rnkPartnerB_", cellCusterBName)
    )
    tmp <- rowSums(result[, interestColumn] < filterCutoff)
    result[tmp != 0, "SCORE"] <- NA
    result <- result[!is.na(result$SCORE), ]
    result <- result[order(result$SCORE, decreasing = TRUE), ]
    result <- data.frame(
      interationPair = rownames(result),
      result,
      stringsAsFactors = FALSE
    )
    columnToremove <- c(
      "partnerA1",
      "partnerA2",
      "partnerA3",
      "partnerB1",
      "partnerB2",
      "partnerB3"
    )
    result <- result[, setdiff(colnames(result), columnToremove)]
    
    if (nrow(result) == 0) {
      stop("No interaction pair exceed the score filterCutoff")
    }
    
    scTHIresult <- list(result, expMat, cellCusterA, cellCusterB)
    names(scTHIresult) <-
      c("result", "expMat", cellCusterAName, cellCusterBName)
    
    ###################### Permutation #########
    if (PValue == TRUE) {
      message("Computing permutation....")
      interaction_table <- interaction_table[rownames(result), ]
      columnToadd <- c(
        "partnerA1",
        "partnerA2",
        "partnerA3",
        "partnerB1",
        "partnerB2",
        "partnerB3"
      )
      result <- data.frame(interaction_table[, columnToadd], result,
                           stringsAsFactors = FALSE
      )
      
      if (length(grep("simple", result$interactionType)) != 0) {
        sample_pairInteraction <-
          rownames(result)[result$interactionType == "simple"][1]
        
        param <- SnowParam(workers = ncore, type = "SOCK")
        if (autocrineEffect == FALSE) {
          paracrineFun <- function(dummy, rankingValues) {
            ddataA <- expMat[, cellCusterA]
            ddataA <-
              apply(ddataA, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataA) <- rownames(expMat)
            ddataA_ <- rankingValues(ddataA)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerA1", "partnerA2",
                "partnerA3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataA_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            ggenes2 <-
              result[sample_pairInteraction, c(
                "partnerB1", "partnerB2",
                "partnerB3"
              )]
            ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
            tmp_genes2 <- ddataA_[ggenes2, , drop = FALSE]
            tmp2 <- colSums(tmp_genes2 == FALSE)
            tmp <- tmp1 + tmp2
            rnk_permuA <- round(sum(tmp == (
              nrow(tmp_genes1) +
                nrow(tmp_genes2)
            )) / ncol(tmp_genes1), digits = 2)
            
            
            ddataB <- expMat[, cellCusterB]
            ddataB <-
              apply(ddataB, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataB) <- rownames(expMat)
            ddataB_ <- rankingValues(ddataB)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerB1", "partnerB2",
                "partnerB3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataB_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            
            ggenes2 <-
              result[sample_pairInteraction, c(
                "partnerA1", "partnerA2",
                "partnerA3"
              )]
            ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
            tmp_genes2 <- ddataB_[ggenes2, , drop = FALSE]
            tmp2 <- colSums(tmp_genes2 == FALSE)
            tmp <- tmp1 + tmp2
            
            rnk_permuB <- round(sum(tmp == (
              nrow(tmp_genes1) +
                nrow(tmp_genes2)
            )) / ncol(tmp_genes1), digits = 2)
            rnkInteraction <- mean(rnk_permuA, rnk_permuB)
            return(rnkInteraction)
          }
          ans <- bplapply(1:nPermu,
                          FUN = paracrineFun,
                          rankingValues = rankingValues,
                          BPPARAM = param
                          
          )
          ans <- unlist(ans)
          permutatedPvalue_simple <- ans
        }
        if (autocrineEffect == TRUE) {
          autocrineFun <- function(dummy, rankingValues) {
            ddataA <- expMat[, cellCusterA]
            ddataA <-
              apply(ddataA, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataA) <- rownames(expMat)
            ddataA_ <- rankingValues(ddataA)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerA1", "partnerA2",
                "partnerA3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataA_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            tmp <- tmp1
            rnk_permuA <-
              round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
                    digits = 2
              )
            
            
            ddataB <- expMat[, cellCusterB]
            ddataB <-
              apply(ddataB, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataB) <- rownames(expMat)
            ddataB_ <- rankingValues(ddataB)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerB1", "partnerB2",
                "partnerB3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataB_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            tmp <- tmp1
            
            rnk_permuB <-
              round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
                    digits = 2
              )
            rnkInteraction <- mean(rnk_permuA, rnk_permuB)
            return(rnkInteraction)
          }
          ans <- bplapply(1:nPermu,
                          FUN = autocrineFun,
                          rankingValues = rankingValues,
                          BPPARAM = param
                          
          )
          ans <- unlist(ans)
          permutatedPvalue_simple <- ans
        }
      }
      
      if (length(grep("complex", result$interactionType)) != 0) {
        sample_pairInteraction <-
          rownames(result)[result$interactionType == "complex"][1]
        
        
        param <- SnowParam(workers = ncore, type = "SOCK")
        if (autocrineEffect == FALSE) {
          paracrineFun2 <- function(dummy, rankingValues) {
            ddataA <- expMat[, cellCusterA]
            ddataA <-
              apply(ddataA, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataA) <- rownames(expMat)
            ddataA_ <- rankingValues(ddataA)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerA1", "partnerA2",
                "partnerA3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataA_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            
            ggenes2 <-
              result[sample_pairInteraction, c(
                "partnerB1", "partnerB2",
                "partnerB3"
              )]
            ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
            tmp_genes2 <- ddataA_[ggenes2, , drop = FALSE]
            tmp2 <- colSums(tmp_genes2 == FALSE)
            
            tmp <- tmp1 + tmp2
            rnk_permuA <- round(sum(tmp == (
              nrow(tmp_genes1) +
                nrow(tmp_genes2)
            )) / ncol(tmp_genes1),
            digits = 2
            )
            
            
            ddataB <- expMat[, cellCusterB]
            ddataB <-
              apply(ddataB, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataB) <- rownames(expMat)
            ddataB_ <- rankingValues(ddataB)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerB1", "partnerB2",
                "partnerB3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataB_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            
            ggenes2 <-
              result[sample_pairInteraction, c(
                "partnerA1", "partnerA2",
                "partnerA3"
              )]
            ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
            tmp_genes2 <- ddataB_[ggenes2, , drop = FALSE]
            tmp2 <- colSums(tmp_genes2 == FALSE)
            tmp <- tmp1 + tmp2
            
            rnk_permuB <- round(sum(tmp == (
              nrow(tmp_genes1) +
                nrow(tmp_genes2)
            )) / ncol(tmp_genes1), digits = 2)
            rnkInteraction <- mean(rnk_permuA, rnk_permuB)
            return(rnkInteraction)
          }
          ans <- bplapply(1:nPermu,
                          FUN = paracrineFun2,
                          rankingValues = rankingValues,
                          BPPARAM = param
                          
          )
          ans <- unlist(ans)
          permutatedPvalue_complex <- ans
        }
        if (autocrineEffect == TRUE) {
          autocrineFun2 <- function(dummy, rankingValues) {
            ddataA <- expMat[, cellCusterA]
            ddataA <-
              apply(ddataA, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataA) <- rownames(expMat)
            ddataA_ <- rankingValues(ddataA)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerA1", "partnerA2",
                "partnerA3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataA_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            tmp <- tmp1
            rnk_permuA <-
              round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
                    digits = 2
              )
            
            
            ddataB <- expMat[, cellCusterB]
            ddataB <-
              apply(ddataB, 2, function(x) {
                sample(x, replace = FALSE)
              })
            rownames(ddataB) <- rownames(expMat)
            ddataB_ <- rankingValues(ddataB)
            
            ggenes1 <-
              result[sample_pairInteraction, c(
                "partnerB1", "partnerB2",
                "partnerB3"
              )]
            ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
            tmp_genes1 <- ddataB_[ggenes1, , drop = FALSE]
            tmp1 <- colSums(tmp_genes1 == TRUE)
            tmp <- tmp1
            
            rnk_permuB <-
              round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
                    digits = 2
              )
            rnkInteraction <- mean(rnk_permuA, rnk_permuB)
            return(rnkInteraction)
          }
          ans <- bplapply(1:nPermu,
                          FUN = autocrineFun2,
                          rankingValues = rankingValues,
                          BPPARAM = param
                          
          )
          ans <- unlist(ans)
          permutatedPvalue_complex <- ans
        }
      }
      
      
      
      r <- apply(result, 1, function (x) {
        if (x["interactionType"] == "simple") {
          ans <-    sum(permutatedPvalue_simple > x["SCORE"]) / nPermu
        }
        if (x["interactionType"] == "complex") {
          ans <- sum(permutatedPvalue_complex > x["SCORE"]) / nPermu
        }
        ans
      })
      SCOREpValue <- unlist(r)
      FDR <- round(p.adjust(SCOREpValue, method = "fdr"), digits = 3)
      result <- data.frame(
        result,
        pValue = SCOREpValue,
        FDR = FDR,
        stringsAsFactors = FALSE
      )
      
      result <- result[result$pValue <= pvalueCutoff, ]
      columnToremove <- c(
        "partnerA1",
        "partnerA2",
        "partnerA3",
        "partnerB1",
        "partnerB2",
        "partnerB3"
      )
      result <- result[, setdiff(colnames(result), columnToremove)]
      scTHIresult$result <- result
      # save(scTHIresult, file = paste0(fileNameBase,
      # "_", cellCusterAName,
      #  "&", cellCusterBName, ".RData"))
    }
    message("Interaction pairs detected:")
    print(rownames(scTHIresult$result))
    return(scTHIresult)
  }
