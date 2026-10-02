#' Replace Dimension Names with Biological Labels
#'
#' Replaces generic model dimension names with biological labels.
#'
#' @keywords internal
switchNames <- function(inResults, infectionLabels, partitionLabels) {
  dimnames(inResults$M_generation)[2:3] <- list(partitionLabels, infectionLabels)
  dimnames(inResults$sigma2_generation)[2:3] <- list(partitionLabels, infectionLabels)
  dimnames(inResults$probs_generation) <- list(partitionLabels, infectionLabels)
  dimnames(inResults$weights_generation)[2:3] <- list(partitionLabels, infectionLabels)
  inResults
}

#' Select by thresholding
#'
#' @keywords internal
select_by_threshold <- function(x, cutoff, f.labels) {
  z <- if (is.numeric(cutoff)) which(x >= cutoff) else which(x == max(x))
  z <- f.labels[z]
  if (length(z) != 1) z <- "NotAssigned"
  z
}

extract_test_components <- function(testData) {
  # Matrix::t keeps sparse layers sparse; consumers densify after gene subsetting
  Y_test <- tryCatch(Matrix::t(testData[["RNA"]]$data.Test), error = function(e) NULL)
  if (is.null(Y_test) && inherits(testData, "Seurat")) {
    Y_test <- tryCatch(
      Matrix::t(SeuratObject::LayerData(testData, assay = "RNA", layer = "data.Test")),
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
      if (is.null(c_obs) && "C.preLabel" %in% colnames(meta)) c_obs <- meta$C.preLabel
      if (is.null(v_obs) && "V.preLabel" %in% colnames(meta)) v_obs <- meta$V.preLabel
    }
  }

  if (is.null(c_obs) || is.null(v_obs)) {
    stop("Missing C.preLabel and/or V.preLabel in testData metadata.", call. = FALSE)
  }
  if (length(c_obs) != nrow(Y_test) || length(v_obs) != nrow(Y_test)) {
    stop("Length mismatch: labels C.preLabel/V.preLabel must match number of test cells.", call. = FALSE)
  }

  list(Y_test = Y_test, c_obs = c_obs, v_obs = v_obs)
}

extract_model_components <- function(mod_results) {
  M <- mod_results$M_generation
  sigma2 <- mod_results$sigma2_generation
  probs <- mod_results$probs_generation

  if (is.null(M) || is.null(sigma2) || is.null(probs)) {
    M <- mod_results$M_list
    sigma2 <- mod_results$sigma2_list
    probs <- mod_results$probs_list
  }
  if (is.null(M) || is.null(sigma2) || is.null(probs)) {
    stop(
      paste0(
        "mod_results is missing model components. Expected either ",
        "M_generation/sigma2_generation/probs_generation or ",
        "M_list/sigma2_list/probs_list."
      ),
      call. = FALSE
    )
  }

  if (is.list(M) || is.list(sigma2) || is.list(probs)) {
    comp_lengths <- c(length(M), length(sigma2), length(probs))
    if (all(comp_lengths == 1L)) {
      M <- M[[1]]
      sigma2 <- sigma2[[1]]
      probs <- probs[[1]]
    } else {
      stop(
        paste0(
          "mod_results contains multiple fitted lambdas. ",
          "Please select one lambda-specific model before calling deTest/getTestWeights."
        ),
        call. = FALSE
      )
    }
  }

  if (is.null(dim(M)) || length(dim(M)) != 3L) {
    stop("Model means must be a 3D array with dims genes x partitions x infection states.", call. = FALSE)
  }
  if (is.null(dim(sigma2)) || length(dim(sigma2)) != 3L) {
    stop("Model variances must be a 3D array with dims genes x partitions x infection states.", call. = FALSE)
  }
  if (is.null(dim(probs)) || length(dim(probs)) != 2L) {
    stop("Model probs must be a 2D matrix with dims partitions x infection states.", call. = FALSE)
  }

  list(M = M, sigma2 = sigma2, probs = probs)
}
