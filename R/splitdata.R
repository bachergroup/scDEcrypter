
#' Split Data into Training and Test Sets by Variant
#'
#' Splits the observed dataset into training ("generation") and test
#' subsets, stratified by variant labels. Samples with missing variant
#' labels are also split proportionally.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param C.obs Integer vector of length \code{n}, giving the observed
#'   condition labels for each sample (may contain \code{NA}).
#' @param V.obs Integer vector of length \code{n}, giving the observed
#'   infection status labels for each sample (may contain \code{NA}).
#' @param seed Integer. Random seed for reproducibility.
#' @param train_frac Numeric in \code{(0,1)}. Proportion of samples
#'   within each variant (and missing group) assigned to the training
#'   set. Defaults to \code{0.5}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{Y_generation}}{Training subset of \code{Y}.}
#'   \item{\code{C.obs_generation}}{Training subset of \code{C.obs}.}
#'   \item{\code{V.obs_generation}}{Training subset of \code{V.obs}.}
#'   \item{\code{Y_test}}{Test subset of \code{Y}.}
#'   \item{\code{C.obs_test}}{Test subset of \code{C.obs}.}
#'   \item{\code{V.obs_test}}{Test subset of \code{V.obs}.}
#'   \item{\code{generation_idx}}{Row indices assigned to training.}
#'   \item{\code{test_idx}}{Row indices assigned to test.}
#' }
#'
#'  @export
splitdata <- function(Y, C.obs, V.obs, seed, train_frac = 0.5) {
  set.seed(seed)
  n <- nrow(Y)
  V.levels <- sort(unique(na.omit(V.obs)))

  generation_idx <- c()
  test_idx <- c()

  for (v in V.levels) {
    idx_v <- which(V.obs == v)
    n_v <- length(idx_v)
    n_train <- round(train_frac * n_v)
    train_idx_v <- sample(idx_v, n_train)
    test_idx_v <- setdiff(idx_v, train_idx_v)

    generation_idx <- c(generation_idx, train_idx_v)
    test_idx <- c(test_idx, test_idx_v)
  }

  idx_na <- which(is.na(V.obs))
  if (length(idx_na) > 0) {
    n_na <- length(idx_na)
    n_train_na <- round(train_frac * n_na)
    train_idx_na <- sample(idx_na, n_train_na)
    test_idx_na <- setdiff(idx_na, train_idx_na)

    generation_idx <- c(generation_idx, train_idx_na)
    test_idx <- c(test_idx, test_idx_na)
  }

  list(
    Y_generation = Y[generation_idx, , drop = FALSE],
    C.obs_generation = C.obs[generation_idx],
    V.obs_generation = V.obs[generation_idx],
    Y_test = Y[test_idx, , drop = FALSE],
    C.obs_test = C.obs[test_idx],
    V.obs_test = V.obs[test_idx],
    generation_idx = generation_idx,
    test_idx = test_idx
  )
}

preprocess_scDEcrypter <- function(seurat_obj, C_obs, V_obs, seed, vs_method = "shifted_log")
{
  vs_method <- match.arg(vs_method)
  set.seed(seed)
  
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("Seurat object must contain an RNA assay.")
  }
  Y <- seurat_obj@assays$RNA$counts
  splitted <- splitdata(t(Y), C_obs, V_obs, seed)
  library(transformGamPoi)
  
  vs_train <- transformGamPoi(
    as.matrix(t(splitted$Y_generation),
              transformation = vs_method,
              size_factors = "poscounts"
    ))
  
  vs_test <- transformGamPoi(
    as.matrix(t(splitted$Y_test),
              transformation = vs_method,
              size_factors = "poscounts"
    ))
  
  Y_gen_stbl <- t(vs_train)
  Y_test_stbl <- t(vs_test)
  rm(vs_train)
  rm(vs_test)
  
  return(list(
    Y_generation = Y_gen_stbl,
    Y_test = Y_test_stbl,
    splitted = splitted
  ))
}
