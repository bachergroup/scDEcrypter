
#' Compute Weights Matrix on Test Data
#'
#' Computes the probability (weight) matrix for test data using the fitted model.
#'
#' @param mod_results A list containing model fitting results from \code{\link{fit_scDEcrypter}}.
#' @param testData A Seurat object or list containing test data.
#' @return The W weight matrix (cells x cell types x viral types).
#' @export
getTestWeights <- function(mod_results, testData) {
  test_components <- extract_test_components(testData)
  Y_test <- test_components$Y_test
  c_obs <- test_components$c_obs
  v_obs <- test_components$v_obs
    M <- mod_results$M_generation
    sigma2 <- mod_results$sigma2_generation
    probs <- mod_results$probs_generation
    
    topGenes <- rownames(M)
    W.test <- E_step(Y_test[, topGenes], c_obs, v_obs, M, probs, sigma2)
    
    mod_results$weights_test <- W.test
    return(mod_results)
}

#' Perform Differential Expression Analysis on Test Data
#'
#' Conducts likelihood ratio testing for differential expression by comparing
#' alternative and null models fitted to test data. Computes log-likelihoods,
#' test statistics, p-values, and log-fold changes for genes across cell type
#' and viral type combinations.
#'
#' @param mod_results A list containing model fitting results from \code{\link{fit_scDEcrypter}}.
#' @param testData A Seurat object or list containing test data with:
#'   \describe{
#'     \item{\code{RNA}}{RNA assay with \code{data.Test} slot containing test counts}
#'     \item{\code{C.preLabel}}{Observed cell type labels for test cells in metadata}
#'     \item{\code{V.preLabel}}{Observed viral status labels for test cells in metadata}
#'   }
#' @param testingGenes Character vector of gene names for which to perform differential
#'   expression testing. Should be a larger set of genes for testing.
#' @param compGroups vector specifying pairwise comparisons to perform. Structure should
#'   indicate which /viral type pairs to compare (e.g., infected vs. uninfected).
#' @param method Character scalar specifying the differential-expression test to use.
#'   Options are \code{"permutation"} (default), \code{"perm"} (alias for
#'   permutation), and \code{"chisq"}.
#' @param nPerm Integer number of within-cell-type permutations to use when
#'   \code{deTest = "permutation"}. Default is \code{1000}.
#' @param mc.cores Integer number of processes to use for permutation testing.
#'   Default is \code{1}.
#'
#' @return A list containing differential expression results:
#'   \describe{
#'     \item{\code{ll.alternative}}{Log-likelihood under alternative (full) model (genes x cell types)}
#'     \item{\code{ll.null}}{Log-likelihood under null model (genes x cell types)}
#'     \item{\code{lrt.stat}}{Likelihood ratio test statistics (genes x cell types)}
#'     \item{\code{pval}}{Permutation or chi-squared p-values (genes x cell types)}
#'     \item{\code{adj.pval}}{FDR-adjusted p-values (genes x cell types)}
#'     \item{\code{logfc}}{Log-fold changes between compared groups}
#'     \item{\code{method}}{The testing method used to compute p-values}
#'   }
#'
#' @details
#' The function fits two models to test data: (1) an alternative model allowing
#' separate means for each comparison group, and (2) a null model with constrained
#' means. A likelihood ratio test compares these nested models for each gene and
#' cell type. By default, p-values are estimated empirically by permuting rows of
#' the weight array within observed cell types, recomputing the likelihood ratio
#' statistic, and comparing permuted statistics to the observed statistic. A
#' chi-squared approximation with 1 degree of freedom is also available and both
#' p-value types are adjusted for multiple testing using the false discovery rate.
#' @importFrom stats pchisq p.adjust
#' @export
deTest <- function(mod_results, testData, 
                testingGenes, compGroups,
                method = c("permutation", "perm", "chisq"),
                nPerm = 1000,
                mc.cores = 1) {
  test.method <- match.arg(method)
  if (identical(test.method, "perm")) {
    test.method <- "permutation"
  }
    
    test_components <- extract_test_components(testData)
    Y_test <- test_components$Y_test
    c_obs <- test_components$c_obs
    v_obs <- test_components$v_obs
    
    if (is.null(mod_results$weights_test)) {
      message("Calculating weights on test set...")
      W.test <- getTestWeights(mod_results, testData)$weights_test
    } else {
      W.test <- mod_results$weights_test
    }
    
    Y_test <- Y_test[, testingGenes]

    message("Calculating observed ...")
    observed_fit <- de_lrt_from_weights(Y_test, W.test, compGroups)
    l1.alt <- observed_fit$ll.alternative
    l0.null <- observed_fit$ll.null
    lrt <- observed_fit$lrt.stat

    if (test.method == "permutation") {
       message("Calculating permutation test statistics ...")
      group_indices <- Filter(length, split(seq_along(c_obs), c_obs))
      permute_apply <- if (nPerm > 1000) parallel::mclapply else pbmclapply
      permuted_lrt <- permute_apply(
        X = seq_len(nPerm),
        FUN = function(iteration) {
          W_perm <- permute_weights_within_groups(W.test, group_indices)
          de_lrt_exceeds_from_weights(Y_test, W_perm, compGroups, lrt)
        },
        mc.cores = mc.cores
      )

      exceed_counts <- matrix(0L,
                              nrow = nrow(lrt),
                              ncol = ncol(lrt),
                              dimnames = dimnames(lrt))
      for (perm_lrt in permuted_lrt) {
        exceed_counts <- exceed_counts + perm_lrt
      }
      pval <- (exceed_counts + 1) / (nPerm + 1)
    } else {
      pval <- apply(lrt, 2, function(x) pchisq(x, df = 1, lower.tail = FALSE))
    }

    pval.adjust <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
    logFC <- compute_logFC(observed_fit$M.alt, compGroups)
    return(list(
        ll.alternative = l1.alt,
        ll.null = l0.null,
        lrt.stat = lrt,
        pval = pval,
        adj.pval = pval.adjust,
        logfc = logFC,
        test.method = test.method
    ))
}

de_lrt_from_weights <- function(Y, W, compGroups) {
    M.alt <- DE_mu(Y = Y, W = W, compStatus = compGroups)
    sigma2.alt <- DE_sigma2(Y, W, M.alt)
    ll.alternative <- approx_complete_data_loglik_fast(Y, M.alt, W, sigma2.alt)

    M.null <- DE_mu_null(Y, W, compGroups)
    sigma2.null <- DE_sigma2(Y, W, M.null)
    ll.null <- approx_complete_data_loglik_fast(Y, M.null, W, sigma2.null)

    list(
      M.alt = M.alt,
      ll.alternative = ll.alternative,
      ll.null = ll.null,
      lrt.stat = -2 * (ll.null - ll.alternative)
    )
}

permute_weights_within_groups <- function(W, group_indices) {
    permuted_index <- seq_len(dim(W)[1])
    for (group_positions in group_indices) {
      if (length(group_positions) > 1) {
        permuted_index[group_positions] <- group_positions[sample.int(length(group_positions))]
      }
    }

    W[permuted_index, , , drop = FALSE]
}



#' Score Cell Type and Viral Status Labels via Thresholding on Test Data
#'
#' Assigns cell type and viral infection status labels to test cells
#' based on conditional weights. Assignments can be made using independent or 
#' joint thresholding approaches.
#'
#' @param mod_results A list containing model fitting results from \code{\link{fit_scDEcrypter}}.
#' @param testData A Seurat object or list containing test data with:
#'   \describe{
#'     \item{\code{RNA}}{RNA assay with \code{data.Test} slot containing test counts}
#'     \item{\code{C.preLabel}}{Observed cell type labels for test cells in metadata}
#'     \item{\code{V.preLabel}}{Observed viral status labels for test cells in metadata}
#'   }
#' @param cutoffInfection Numeric threshold in \code{[0,1]} for viral status prediction,
#'   or \code{"max"} to select the infection status with highest weight. Default: \code{0.9}.
#' @param cutoffPartition Numeric threshold in \code{[0,1]} for cell type prediction,
#'   or \code{"max"} to select the cell type with highest weight. Default: \code{"max"}.
#' @param independentThres Logical. If \code{TRUE} (default), cell type and viral status
#'   are predicted independently by marginalizing over the other dimension. If \code{FALSE},
#'   joint predictions are made using the maximum weight across both dimensions.
#'
#' @return The input \code{testData} object with two additional columns added to metadata:
#'   \describe{
#'     \item{\code{Infection_pred}}{Predicted viral infection status labels (Unknown = uninfected/uncertain)}
#'     \item{\code{Partition_pred}}{Predicted partition variable (e.g., cell type) labels (Unknown = unclassified/uncertain)}
#'   }
#'   A value of Unknown indicates the cell could not be confidently assigned (either multiple
#'   dimensions exceeded threshold or none exceeded threshold).
#'
#' @details
#' When \code{independentThres = TRUE}, the function sums weights across the alternative
#' dimension (e.g., summing across viral types when predicting cell type). This marginalizes
#' out uncertainty in the other dimension.
#'
#' When \code{independentThres = FALSE}, the maximum weight is taken across each dimension
#' separately, then thresholding is applied. This provides joint predictions.
#'
#' Numeric thresholds allow for conservative (higher threshold) or liberal (lower threshold)
#' predictions. The special value \code{"max"} assigns each cell to the category with
#' highest weight, which always produces a definitive assignment unless there are ties (returns as Unknown).
#'
#' @export
thresholdScoring <- function(mod_results, testData, cutoffInfection = .9, 
                  cutoffPartition = "max", independentThres = TRUE) {
    
    Y_test <- t(testData[["RNA"]]$data.Test)
    
    if (is.null(mod_results$weights_test)) {
      message("Calculating weights on test set...")
      W.test <- getTestWeights(mod_results, testData)$weights_test
    } else {
      W.test <- mod_results$weights_test
    }
    
    c.labs = dimnames(W.test)[[2]]
    v.labs = dimnames(W.test)[[3]]
    message("Predicting on test set...")
    
    ## Based on user-input (sum is indep., max is not)
    reduce_fn <- if (independentThres) sum else max
    
    W_byC <- apply(W.test, c(1, 2), reduce_fn, na.rm = TRUE)
    Partition_pred <- apply(W_byC, 1, function(x) select_by_threshold(x, cutoffPartition, c.labs))

    W_byV <- apply(W.test, c(1, 3), reduce_fn, na.rm = TRUE)
    Infection_pred <- unlist(apply(W_byV, 1, function(x) select_by_threshold(x, cutoffInfection, v.labs)))
    
    testData$pred_InfectionState <- Infection_pred
    testData$pred_PartitionState <- Partition_pred
    
    return(testData)
}


#' Replace Dimension Names with Biological Labels
#'
#' Replaces the generic dimension names (e.g., "partitionStatus_1", "infectionStatus_2")
#' in model results with meaningful biological labels for improved interpretability
#' and downstream analysis.
#'
#' @param inResults A list containing model fitting results from \code{\link{fit_scDEcrypter}}.
#' @param infectionLabels Character vector of length equal to the number of viral types,
#'   specifying biological names for each infection/variant status. 
#'   Example: \code{c("Uninfected", "Infected", "Bystander", etc)}
#' @param partitionLabels Character vector of length equal to the number of cell types,
#'   specifying biological names for each cell type cluster.
#'   Example: \code{c("T_cells", "B_cells", "Monocytes", etc)}
#'
#' @return The input \code{inResults} list with updated dimension names:
#'   \describe{
#'     \item{\code{M_generation}}{Dimension names: genes x \code{partitionLabels} x \code{infectionLabels}}
#'     \item{\code{sigma2_generation}}{Dimension names: genes x \code{partitionLabels} x \code{infectionLabels}}
#'     \item{\code{probs_generation}}{Dimension names: \code{partitionLabels} x \code{infectionLabels}}
#'     \item{\code{weights_generation}}{Dimension names: cells x \code{partitionLabels} x \code{infectionLabels}}
#'   }
#'
#' @details
#' This function is useful for converting generic cluster labels generated during
#' model fitting into interpretable biological labels after validation and annotation.
#' The updated dimension names will be preserved through downstream analyses and
#' improve readability of results tables and visualizations.
#'
#' @keywords internal
switchNames <- function(inResults, infectionLabels, partitionLabels) {
    
    dimnames(inResults$M_generation)[2:3] <- list(partitionLabels, infectionLabels)
    dimnames(inResults$sigma2_generation)[2:3]  <- list(partitionLabels, infectionLabels)
    dimnames(inResults$probs_generation) <- list(partitionLabels, infectionLabels)
    dimnames(inResults$weights_generation)[2:3] <- list(partitionLabels, infectionLabels)
    
    return(inResults)
}


#' Select by thresholding
#'
#' @keywords internal
select_by_threshold <- function(x, cutoff, f.labels) {
    z <- if (is.numeric(cutoff)) which(x >= cutoff) else which(x == max(x))
    z <- f.labels[z]
    if (length(z) != 1) z <- "NotAssigned"
    return(z)
}


extract_test_components <- function(testData) {
    Y_test <- tryCatch(t(testData[["RNA"]]$data.Test), error = function(e) NULL)
    if (is.null(Y_test) && inherits(testData, "Seurat")) {
      Y_test <- tryCatch(
        t(SeuratObject::LayerData(testData, assay = "RNA", layer = "data.Test")),
        error = function(e) NULL
      )
    }
    if (is.null(Y_test)) {
      stop(
        "Unable to read test expression matrix. Expected RNA layer/slot 'data.Test'.",
        call. = FALSE
      )
    }

    c_obs <- tryCatch(testData$C.preLabel, error = function(e) NULL)
    v_obs <- tryCatch(testData$V.preLabel, error = function(e) NULL)

    if ((is.null(c_obs) || is.null(v_obs)) && inherits(testData, "Seurat")) {
      meta <- tryCatch(testData@meta.data, error = function(e) NULL)
      if (!is.null(meta)) {
        if (is.null(c_obs) && "C.preLabel" %in% colnames(meta)) {
          c_obs <- meta$C.preLabel
        }
        if (is.null(v_obs) && "V.preLabel" %in% colnames(meta)) {
          v_obs <- meta$V.preLabel
        }
      }
    }

    if (is.null(c_obs) || is.null(v_obs)) {
      stop(
        "Missing C.preLabel and/or V.preLabel in testData metadata.",
        call. = FALSE
      )
    }

    if (length(c_obs) != nrow(Y_test) || length(v_obs) != nrow(Y_test)) {
      stop(
        "Length mismatch: labels C.preLabel/V.preLabel must match number of test cells.",
        call. = FALSE
      )
    }

    list(Y_test = Y_test, c_obs = c_obs, v_obs = v_obs)
}