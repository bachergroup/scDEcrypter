
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
    model_components <- extract_model_components(mod_results)
    M <- model_components$M
    sigma2 <- model_components$sigma2
    probs <- model_components$probs
    
    topGenes <- rownames(M)
    W.test <- E_step(Y_test[, topGenes, drop = FALSE], c_obs, v_obs, M, probs, sigma2)
    
    mod_results$weights_test <- W.test
    return(mod_results)
}

#' Score Cell Type and Viral Status Labels via Thresholding on Test Data
#'
#' Assigns cell type and viral infection status labels to test cells based on
#' conditional weights.
#'
#' @param mod_results A list containing model fitting results from \code{\link{fit_scDEcrypter}}.
#' @param testData A Seurat object or list containing test data.
#' @param cutoffInfection Numeric threshold in \code{[0,1]}, or \code{"max"}.
#' @param cutoffPartition Numeric threshold in \code{[0,1]}, or \code{"max"}.
#' @param independentThres If \code{TRUE}, marginalize over the other dimension.
#' @return The input test data with predicted infection and partition states.
#' @export
thresholdScoring <- function(mod_results, testData, cutoffInfection = .9,
                             cutoffPartition = "max", independentThres = TRUE) {
  if (is.null(mod_results$weights_test)) {
    message("Calculating weights on test set...")
    W.test <- getTestWeights(mod_results, testData)$weights_test
  } else {
    W.test <- mod_results$weights_test
  }

  c.labs <- dimnames(W.test)[[2]]
  v.labs <- dimnames(W.test)[[3]]
  reduce_fn <- if (independentThres) sum else max

  W_byC <- apply(W.test, c(1, 2), reduce_fn, na.rm = TRUE)
  Partition_pred <- apply(W_byC, 1, function(x) {
    select_by_threshold(x, cutoffPartition, c.labs)
  })
  W_byV <- apply(W.test, c(1, 3), reduce_fn, na.rm = TRUE)
  Infection_pred <- unlist(apply(W_byV, 1, function(x) {
    select_by_threshold(x, cutoffInfection, v.labs)
  }))

  testData$pred_InfectionState <- Infection_pred
  testData$pred_PartitionState <- Partition_pred
  testData
}

