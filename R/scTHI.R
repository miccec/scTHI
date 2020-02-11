scTHI.score <-
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
    
        #' scTHI.score
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
        #' ####################### example of scTHI.score
        #' data(scExample)
        #' result <-  scTHI.score(scExample,
        #'       cellCusterA = colnames(scExample)[1:30],
        #'       cellCusterB = colnames(scExample)[31:100],
        #'       cellCusterAName = "ClusterA",
        #'       cellCusterBName = "ClusterB", filterCutoff = 0,
        #'      pvalueCutoff = 1, nPermu = 100, ncore = 8)
        #'
        #' @return A list of results, with four items: result (data.frame),
        #'   expMat (matrix), clusterA (character),  clusterA (character)
        #' @export
        #' scTHI.score
        ############### check all interaction genes in expMat ################
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

        #################### scTHI score ##############################
        message(paste(
            "Computing score for",
            nrow(interaction_table),
            "interaction pairs..."
        ))

        ddataA <- expMat[, cellCusterA]
        n <- round(nrow(ddataA) * topRank / 100)
        ddataA_ <- apply(ddataA, 2, function(x) {
            rank(-x)
        })
        ddataA_ <- apply(ddataA_, 2, function(x) {
            x <= n
        })
        rnkA <- rep(0, nrow(interaction_table))
        expValueA <- rep(0, nrow(interaction_table))
        names(rnkA) <- rownames(interaction_table)
        names(expValueA) <- rownames(interaction_table)

        if (autocrineEffect == FALSE) {
            for (i in 1:nrow(interaction_table)) {
                ggenes1 <-
                    interaction_table[i, c("partnerA1", "partnerA2", 
                                           "partnerA3")]
                ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
                tmp_genes1 <- ddataA_[ggenes1, , drop = FALSE]
                tmp1 <- colSums(tmp_genes1 == TRUE)

                ggenes2 <-
                    interaction_table[i, c("partnerB1", "partnerB2", 
                                           "partnerB3")]
                ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
                tmp_genes2 <- ddataA_[ggenes2, , drop = FALSE]
                tmp2 <- colSums(tmp_genes2 == FALSE)
                tmp <- tmp1 + tmp2

                rnkA[i] <- round(sum(tmp == (nrow(tmp_genes1) +
                    nrow(tmp_genes2))) / ncol(tmp_genes1), digits = 2)
                expValueA[i] <-
                    paste(round(rowMeans(ddataA[ggenes1, cellCusterA,
                        drop = FALSE
                    ]), digits = 2), collapse = " # ")
            }
        }
        if (autocrineEffect == TRUE) {
            for (i in 1:nrow(interaction_table)) {
                ggenes1 <-
                    interaction_table[i, c("partnerA1", "partnerA2", 
                                           "partnerA3")]
                ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
                tmp_genes1 <- ddataA_[ggenes1, , drop = FALSE]
                tmp1 <- colSums(tmp_genes1 == TRUE)
                tmp <- tmp1

                rnkA[i] <-
                    round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
                        digits = 2
                    )
                expValueA[i] <-
                    paste(round(rowMeans(ddataA[ggenes1, cellCusterA,
                        drop = FALSE
                    ]), digits = 2), collapse = " # ")
            }
        }

        resultA <- data.frame(
            rnkPartnerA = rnkA,
            expValueA = expValueA,
            stringsAsFactors = FALSE
        )
        colnames(resultA) <-
            paste0(colnames(resultA), "_", cellCusterAName)
        print(paste("Computed", i, "ranked values for partner A"))


        ddataB <- expMat[, cellCusterB]
        n <- round(nrow(ddataB) * topRank / 100)
        ddataB_ <- apply(ddataB, 2, function(x) {
            rank(-x)
        })
        ddataB_ <- apply(ddataB_, 2, function(x) {
            x <= n
        })
        rnkB <- rep(0, nrow(interaction_table))
        expValueB <- rep(0, nrow(interaction_table))
        names(rnkB) <- rownames(interaction_table)
        names(expValueB) <- rownames(interaction_table)

        if (autocrineEffect == FALSE) {
            for (i in 1:nrow(interaction_table)) {
                ggenes1 <-
                    interaction_table[i, c("partnerB1", "partnerB2", 
                                           "partnerB3")]
                ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
                tmp_genes1 <- ddataB_[ggenes1, , drop = FALSE]
                tmp1 <- colSums(tmp_genes1 == TRUE)

                ggenes2 <-
                    interaction_table[i, c("partnerA1", "partnerA2", 
                                           "partnerA3")]
                ggenes2 <- as.vector(ggenes2[!is.na(ggenes2)])
                tmp_genes2 <- ddataB_[ggenes2, , drop = FALSE]
                tmp2 <- colSums(tmp_genes2 == FALSE)
                tmp <- tmp1 + tmp2

                rnkB[i] <- round(sum(tmp == (nrow(tmp_genes1) +
                    nrow(tmp_genes2))) / ncol(tmp_genes1), digits = 2)
                expValueB[i] <-
                    paste(round(rowMeans(ddataB[ggenes1, cellCusterB,
                        drop = FALSE
                    ]), digits = 2), collapse = " # ")
            }
        }
        if (autocrineEffect == TRUE) {
            for (i in 1:nrow(interaction_table)) {
                ggenes1 <-
                    interaction_table[i, c("partnerB1", "partnerB2", 
                                           "partnerB3")]
                ggenes1 <- as.vector(ggenes1[!is.na(ggenes1)])
                tmp_genes1 <- ddataB_[ggenes1, , drop = FALSE]
                tmp1 <- colSums(tmp_genes1 == TRUE)
                tmp <- tmp1

                rnkB[i] <-
                    round(sum(tmp == nrow(tmp_genes1)) / ncol(tmp_genes1),
                        digits = 2
                    )
                expValueB[i] <-
                    paste(round(rowMeans(ddataB[ggenes1, cellCusterB,
                        drop = FALSE
                    ]), digits = 2), collapse = " # ")
            }
        }

        resultB <- data.frame(
            rnkPartnerB = rnkB,
            expValueB = expValueB,
            stringsAsFactors = FALSE
        )
        colnames(resultB) <-
            paste0(colnames(resultB), "_", cellCusterBName)
        print(paste("Computed", i, "ranked values for partner B"))

        SCORE <- rowMeans(cbind(rnkA, rnkB))
        result <- data.frame(interaction_table, resultA, resultB, SCORE)

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
                    paracrineFun <- function(dummy) {
                        ddataA <- expMat[, cellCusterA]
                        ddataA <-
                            apply(ddataA, 2, function(x) {
                                sample(x, replace = FALSE)
                            })
                        rownames(ddataA) <- rownames(expMat)

                        n <- round(nrow(ddataA) * topRank / 100)
                        ddataA_ <- apply(ddataA, 2, function(x) {
                            rank(-x)
                        })
                        ddataA_ <- apply(ddataA_, 2, function(x) {
                            x <= n
                        })
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
                        n <- round(nrow(ddataB) * topRank / 100)
                        ddataB_ <- apply(ddataB, 2, function(x) {
                            rank(-x)
                        })
                        ddataB_ <- apply(ddataB_, 2, function(x) {
                            x <= n
                        })

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
                        BPPARAM = param
                    )
                    ans <- unlist(ans)
                    permutatedPvalue_simple <- ans
                }
                if (autocrineEffect == TRUE) {
                    autocrineFun <- function(dummy) {
                        ddataA <- expMat[, cellCusterA]
                        ddataA <-
                            apply(ddataA, 2, function(x) {
                                sample(x, replace = FALSE)
                            })
                        rownames(ddataA) <- rownames(expMat)

                        n <- round(nrow(ddataA) * topRank / 100)
                        ddataA_ <- apply(ddataA, 2, function(x) {
                            rank(-x)
                        })
                        ddataA_ <- apply(ddataA_, 2, function(x) {
                            x <= n
                        })

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
                        n <- round(nrow(ddataB) * topRank / 100)
                        ddataB_ <- apply(ddataB, 2, function(x) {
                            rank(-x)
                        })
                        ddataB_ <- apply(ddataB_, 2, function(x) {
                            x <= n
                        })

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
                    paracrineFun2 <- function(dummy) {
                        ddataA <- expMat[, cellCusterA]
                        ddataA <-
                            apply(ddataA, 2, function(x) {
                                sample(x, replace = FALSE)
                            })
                        rownames(ddataA) <- rownames(expMat)

                        n <- round(nrow(ddataA) * topRank / 100)
                        ddataA_ <- apply(ddataA, 2, function(x) {
                            rank(-x)
                        })
                        ddataA_ <- apply(ddataA_, 2, function(x) {
                            x <= n
                        })

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
                        n <- round(nrow(ddataB) * topRank / 100)
                        ddataB_ <- apply(ddataB, 2, function(x) {
                            rank(-x)
                        })
                        ddataB_ <- apply(ddataB_, 2, function(x) {
                            x <= n
                        })

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
                        BPPARAM = param
                    )
                    ans <- unlist(ans)
                    permutatedPvalue_complex <- ans
                }
                if (autocrineEffect == TRUE) {
                    autocrineFun2 <- function(dummy) {
                        ddataA <- expMat[, cellCusterA]
                        ddataA <-
                            apply(ddataA, 2, function(x) {
                                sample(x, replace = FALSE)
                            })
                        rownames(ddataA) <- rownames(expMat)

                        n <- round(nrow(ddataA) * topRank / 100)
                        ddataA_ <- apply(ddataA, 2, function(x) {
                            rank(-x)
                        })
                        ddataA_ <- apply(ddataA_, 2, function(x) {
                            x <= n
                        })

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
                        n <- round(nrow(ddataB) * topRank / 100)
                        ddataB_ <- apply(ddataB, 2, function(x) {
                            rank(-x)
                        })
                        ddataB_ <- apply(ddataB_, 2, function(x) {
                            x <= n
                        })

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
                        BPPARAM = param
                    )
                    ans <- unlist(ans)
                    permutatedPvalue_complex <- ans
                }
            }

            SCOREpValue <- c()
            for (k in 1:nrow(result)) {
                if (result[k, "interactionType"] == "simple") {
                    SCOREpValue <- c(
                        SCOREpValue,
                        sum(permutatedPvalue_simple > result[k, "SCORE"]) / nPermu
                    )
                }
                if (result[k, "interactionType"] == "complex") {
                    SCOREpValue <- c(
                        SCOREpValue,
                        sum(permutatedPvalue_complex > result[k, "SCORE"]) / nPermu
                    )
                }
            }

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


scTHI.plotResult <- function(scTHIresult,
                             cexNames = 0.8,
                             plotType = c("score", "pair"),
                             nRes = NULL) {
    #' scTHI.plotResult
    #'
    #' Creates barplots of scTHI.score results.
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
    #' result <-  scTHI.score(scExample,
    #'                        cellCusterA = colnames(scExample)[1:30],
    #'                        cellCusterB = colnames(scExample)[31:100],
    #'                        cellCusterAName = "ClusterA",
    #'                        cellCusterBName = "ClusterB", filterCutoff = 0,
    #'                        pvalueCutoff = 1, nPermu = 100, ncore = 8)
    #'
    #' scTHI.plotResult(result, plotType = "score")
    #' scTHI.plotResult(result, plotType = "pair")
    #' @return None
    #' @export
    #'
    #' scTHI.plotResult


    result <- scTHIresult$result
    if (!is.null(nRes)) {
        result <- result[1:nRes, ]
    }


    if (plotType == "score") {
        tmp <- as.matrix(result[, c("SCORE", "pValue")])
        pvalues <- tmp[, "pValue"]
        for (i in 1:length(pvalues)) {
            if (tmp[i, "pValue"] > 0.05) {
                pvalues[i] <- "ns"
            }
            if (tmp[i, "pValue"] <= 0.05) {
                pvalues[i] <- "*"
            }
            if (tmp[i, "pValue"] <= 0.01) {
                pvalues[i] <- "**"
            }
            if (tmp[i, "pValue"] <= 0.001) {
                pvalues[i] <- "***"
            }
            if (tmp[i, "pValue"] <= 0.0001) {
                pvalues[i] <- "****"
            }
        }

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


scTHI.runTsne <- function(scTHIresult) {
    #' scTHI.runTsne
    #'
    #' Runs t-SNE dimensionality reduction on selected features based on Rtsne
    #' package.
    #' @param scTHIresult scTHI object.
    #' @examples
    #' data(scExample)
    #' result <-  scTHI.score(scExample,
    #'                        cellCusterA = colnames(scExample)[1:30],
    #'                        cellCusterB = colnames(scExample)[31:100],
    #'                        cellCusterAName = "ClusterA",
    #'                        cellCusterBName = "ClusterB", filterCutoff = 0,
    #'                        pvalueCutoff = 1, nPermu = 100, ncore = 8)
    #' result <- scTHI.runTsne(result)
    #' @return The same object as scTHI.score with a fifth item tsneData
    #'   (data.frame)
    #' @export
    #'
    #' scTHI.runTsne


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
    set.seed(1) ### for reproducibility
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


scTHI.plotCluster <- function(scTHIresult,
                              cexPoint = 0.8,
                              legendPos = c(
                                  "topleft", "topright",
                                  "bottomright", "bottomleft"
                              )) {
    #' scTHI.plotCluster
    #'
    #' Graphs the output of scTHI.runTsne, labeling cells by clusters.
    #' @param scTHIresult scTHI object.
    #' @param cexPoint Set the point size.
    #' @param legendPos Character string to custom the legend position.
    #' @examples
    #' data(scExample)
    #' result <-  scTHI.score(scExample,
    #'                        cellCusterA = colnames(scExample)[1:30],
    #'                        cellCusterB = colnames(scExample)[31:100],
    #'                        cellCusterAName = "ClusterA",
    #'                        cellCusterBName = "ClusterB", filterCutoff = 0,
    #'                        pvalueCutoff = 1, nPermu = 100, ncore = 8)
    #' result <- scTHI.runTsne(result)
    #' scTHI.plotCluster(result)
    #' @return None
    #' @export
    #'
    #' scTHI.plotCluster


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
    for (i in 1:length(Colors)) {
        ClusterColors[scTHIresult[[names(Colors)[i]]]] <- Colors[i]
    }
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


scTHI.plotPairs <-
    function(scTHIresult,
             cexPoint = 0.8,
             interactionToplot) {
        #' scTHI.plotPairs
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
        #' data(scExample)
        #' result <-  scTHI.score(scExample,
        #'                  cellCusterA = colnames(scExample)[1:30],
        #'                  cellCusterB = colnames(scExample)[31:100],
        #'                  cellCusterAName = "ClusterA",
        #'                  cellCusterBName = "ClusterB", filterCutoff = 0,
        #'                  pvalueCutoff = 1, nPermu = 100, ncore = 8)
        #' result <- scTHI.runTsne(result)
        #' scTHI.plotPairs(result,interactionToplot = "CXCL12_CD4")

        #' @return None
        #' @export
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

        list_colori_A <- vector("list", length(genesToplotA))
        names(list_colori_A) <- genesToplotA
        for (i in 1:length(genesToplotA)) {
            tmp_exp <-
                t(expMat[genesToplotA[i], colnames(expMat), drop = FALSE])
            mean_exp_samples <- mean(tmp_exp)

            ##### range gene x gene
            range_exp <- range(tmp_exp)
            bre <- seq(
                min(tmp_exp[, genesToplotA[i]]) - 0.01,
                max(tmp_exp[, genesToplotA[i]]) + 0.01, 0.01
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
            list_colori_A[[i]] <- tmp
        }

        list_colori_B <- vector("list", length(genesToplotB))
        names(list_colori_B) <- genesToplotB
        for (i in 1:length(genesToplotB)) {
            tmp_exp <-
                t(expMat[genesToplotB[i], colnames(expMat), drop = FALSE])
            mean_exp_samples <- mean(tmp_exp)

            ##### range gene x gene
            range_exp <- range(tmp_exp)
            bre <- seq(min(tmp_exp[, genesToplotB[i]]) - 0.01, max(tmp_exp[
                , genesToplotB[i]
            ]) + 0.01, 0.01)
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
            list_colori_B[[i]] <- tmp
        }

        ################### plot
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
        for (i in 1:length(list_colori_A)) {
            colori <- list_colori_A[[i]]$colori
            range_exp <- round(list_colori_A[[i]]$range_exp, digits = 2)
            mean_exp_samples <-
                round(list_colori_A[[i]]$mean_exp_samples, digits = 2)
            plot(
                tsneData[, 1:2],
                pch = 16,
                cex = cexPoint,
                col = colori,
                main = genesToplotA[i],
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
        }

        for (i in 1:length(list_colori_B)) {
            colori <- list_colori_B[[i]]$colori
            range_exp <- round(list_colori_B[[i]]$range_exp, digits = 2)
            mean_exp_samples <-
                round(list_colori_B[[i]]$mean_exp_samples, digits = 2)
            plot(
                tsneData[, 1:2],
                pch = 16,
                cex = cexPoint,
                col = colori,
                main = genesToplotB[i],
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
        }
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
    }



TME.classification <- function(expMat,
                               minLenGeneSet = 10,
                               alternative = "two.sided",
                               pvalFilter = FALSE,
                               fdrFilter = TRUE,
                               pvalCutoff = 0.01,
                               nesCutoff = 0.58,
                               nNES = 1) {
    #' TME.classification
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
    #' data(scExample)
    #' Class <- TME.classification(scExample)
    #' @return A list with two items: Class (character) and ClassLegend
    #'   (character)
    #' @export
    #'
    #' TME.classification


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
    signatures <- signatures[tmp > minLenGeneSet]

    ans <- vector("list", length(signatures))
    names(ans) <- names(signatures)

    for (i in 1:length(signatures)) {
        geneSet <- signatures[[i]]
        gs <- geneSet[which(geneSet %in% rownames(E))]
        outside_gs <- setdiff(rownames(E), gs)
        nx <- length(gs)
        ny <- length(outside_gs)


        Tstat <- colSums(E_[outside_gs, ])
        Ustat <- nx * ny + ny * (ny + 1) / 2 - Tstat
        mu <- nx * ny / 2
        sigma <- sqrt(mu * (nx + ny + 1) / 6)
        zValue <- Ustat - mu
        correction <- sapply(zValue, function(x) {
            switch(
                alternative,
                two.sided = sign(x) * 0.5,
                greater = 0.5,
                less = -0.5
            )
        })

        zValue <- (zValue - correction) / sigma
        pValue <- sapply(zValue, function(x) {
            switch(
                alternative,
                less = 1 - pnorm(x),
                greater = pnorm(-x),
                two.sided = 2 * pnorm(-abs(x))
            )
        })


        nes <- Ustat / nx / ny
        pu <- nes / (1 - nes)
        log.pu <- log2(pu)
        names(Ustat) <- colnames(E_)
        names(pValue) <- colnames(E_)
        names(nes) <- colnames(E_)
        names(pu) <- colnames(E_)
        names(log.pu) <- colnames(E_)
        ans[[i]] <- list(
            statistic = Ustat,
            p.value = pValue,
            nes = nes,
            pu = pu,
            log.pu = log.pu
        )
    }

    NES <- sapply(ans, function(x) {
        cbind(x$log.pu)
    })
    NES <- t(NES)
    colnames(NES) <- colnames(expMat)
    pValue <- sapply(ans, function(x) {
        cbind(x$p.value)
    })
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
    print(sort(table(Class), decreasing = TRUE))
    Classification <- list(Class, ClassLegend)
    names(Classification) <- c("Class", "ClassLegend")
    return(Classification)
}



TME.plot <- function(tsneData, Class, cexPoint = 0.8) {
    #' TME.plot
    #'
    #' Generates a plot on the t-SNE coordinates, labeling cells by TME
    #' classification.
    #'
    #' @param tsneData X and y coordinates of points in the plot.
    #' @param Class Object returned by TME.classification function.
    #' @param cexPoint Set the point size.
    #' @examples
    #' data(scExample)
    #' result <-  scTHI.score(scExample,
    #'            cellCusterA = colnames(scExample)[1:30],
    #'            cellCusterB = colnames(scExample)[31:100],
    #'            cellCusterAName = "ClusterA",
    #'            cellCusterBName = "ClusterB", filterCutoff = 0,
    #'            pvalueCutoff = 1, nPermu = 100, ncore = 8)
    #' result <- scTHI.runTsne(result)
    #' Class <- TME.classification(scExample)
    #' TME.plot(tsneData = result$tsneData, Class)
    #' @return  None
    #' @export
    #'
    #' TME.plot

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
