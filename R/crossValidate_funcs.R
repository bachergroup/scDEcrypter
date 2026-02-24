
#' Cross-validation for Lambda Selection
#'
#' Performs K-fold cross-validation to select the best penalty parameter
#' \eqn{\lambda} for the \code{MultiwayMixture} model. Each fold splits the
#' data into training and validation sets, fits the model on training data
#' across candidate \eqn{\lambda} values, and evaluates log-likelihood on
#' validation data. The best \eqn{\lambda} is chosen by maximizing the
#' average validation log-likelihood across folds.
#'
#' @param Data Numeric matrix of dimension \code{n x p}, training data.
#'   Each row is a sample (cell), each column is a feature (e.g. gene).
#' @param lambda.cands Numeric vector. Candidate \eqn{\lambda} values to be
#'   cross-validated.
#' @param max.iter Integer. Maximum number of EM iterations for
#'   \code{MultiwayMixture}.
#' @param tol Numeric. Convergence tolerance for EM updates.
#' @param c_star Integer. Prespecified number of levels for the partitioning variable (e.g., cell-type).
#' @param v_star Integer. Prespecified number of levels for infection.
#' @param NCORES Integer. Number of cores for parallel computation.
#' @param K Integer. Number of cross-validation folds. Defaults to 5.
#'
#' @importFrom stats na.omit
#' @importFrom caret createFolds
#' @importFrom pbmcapply pbmclapply
#' @return A list with components:
#' \describe{
#'   \item{\code{best_lambda}}{The \eqn{\lambda} value with highest average
#'   validation log-likelihood.}
#'   \item{\code{lambda_log_likelihoods}}{Matrix of log-likelihoods, rows =
#'   \code{lambda.vec}, columns = folds.}
#'   \item{\code{avg_log_likelihoods}}{Vector of mean log-likelihoods across folds.}
#' }
#' @export

cross_validate_lambda <- function(Data, 
                                  lambda.cands = c(1e2, 1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5), 
                                  max.iter = 20, tol = 1e-8, 
                                  c_star, v_star, 
                                  NCORES, K = 5) {
    
    if (length(lambda.cands) < 2) {
        stop("Must enter more than one candidate value for lambda to run cross validation.")
      }
    
    v_obs <- tryCatch(
              Data[["V.preLabel", drop=TRUE]],
              error = function(e) {
                stop("Column 'V.preLabel' not found in Data")
              })
    V.levels <- sort(unique(na.omit(v_obs)))
    folds_list <- vector("list", K)
    
    for (v in V.levels) {
        idx_v <- names(which(v_obs == v))
        fold_indices_v <- caret::createFolds(y=seq_along(idx_v), k = K, 
                                      list = TRUE, returnTrain = FALSE)
        folds_list <- lapply(fold_indices_v, function(x) idx_v[x])
    }
    
    idx_na <- names(which(is.na(v_obs)))
    if (length(idx_na) > 0) {
        fold_indices_na <- caret::createFolds(seq_along(idx_na), k = K, 
                                       list = TRUE, returnTrain = FALSE)
        folds_list <- lapply(seq_len(K), function(x) {
            return(c(folds_list[[x]], idx_na[fold_indices_na[[x]]]))
        })
    }
    
   if (NCORES %% K == 0) {
      main_cores <- K
      sub_cores_per_main <- rep(floor(NCORES / K), K)
   } else {
       main_cores = NCORES
       sub_cores_per_main = 1
   }
  
    lambda_log_likelihoods <- pbmcapply::pbmclapply(seq_len(K), runFOLDS, 
                                       fold_indices = folds_list, 
                                       trainData = Data,
                                       max.iter = max.iter, 
                                       tol = tol, 
                                       c_star = c_star, v_star = v_star, 
                                       lambda.cands = lambda.cands, 
                                       foldCores = sub_cores_per_main,
                                       mc.cores = main_cores)
    
    lambda_log_likelihoods <- do.call(cbind, lambda_log_likelihoods)
    colnames(lambda_log_likelihoods) <- paste0("Fold_", 1:K)
    rownames(lambda_log_likelihoods) <- lambda.cands
    
    avg_log_likelihoods <- rowMeans(lambda_log_likelihoods)
    best_lambda <- lambda.cands[which.max(avg_log_likelihoods)]
    cat("Best lambda selected: ", best_lambda, "\n")
    
    return(list(
        "best_lambda" = best_lambda,
        "lambda_log_likelihoods" = lambda_log_likelihoods,
        "avg_log_likelihoods" = avg_log_likelihoods
    ))
}

#' @keywords internal
fit_single_lambda <- function(lambda, Data, c_star, v_star, max.iter, tol) {

  fit_scDEcrypter(Data = Data, 
                  max.iter = max.iter, tol = tol, 
                  c_star = c_star, v_star = v_star,
                  lambda.vec = lambda)
}

#' @keywords internal
runFOLDS <- function(x, fold_indices, trainData, 
                     max.iter, tol, c_star, v_star,
                     lambda.cands, foldCores) {
    
    message("[Fold ", x, "] Starting...")
    validation_idx <- fold_indices[[x]]
    training_idx <- setdiff(colnames(trainData), validation_idx)
    
    train_trainData <- trainData[,training_idx]
    
    # Perform lambda selection on training set
		results <- pbmcapply::pbmclapply(lambda.cands, 
												  fit_single_lambda,
												  Data = train_trainData,
												  c_star = c_star,
												  v_star = v_star,
												  max.iter = max.iter,
												  tol = tol,
												  mc.cores = foldCores
		)
    
    message("[Fold ", x, "] Fitted model, evaluating lambdas...")
    
    # Evaluate each lambda on the validation set
    Y_val <- t(trainData[["RNA"]]$data.Generation[,validation_idx])
    c_obs_val <- trainData$C.preLabel[validation_idx]
    v_obs_val <- trainData$V.preLabel[validation_idx]
    
    log_likelihoods <- sapply(seq_along(lambda.cands), function(y) {
        observed.data.loglik(Y = Y_val, 
                             M = results[[y]]$M_generation, 
                             sigma2 = results[[y]]$sigma2_generation,
                             probs = results[[y]]$probs_generation, 
                             c_obs = c_obs_val, 
                             v_obs = v_obs_val)
        
    })
   names(log_likelihoods) <- lambda.cands

   message("[Fold ", x, "] Complete")
   return(log_likelihoods)
}



#' Observed-Data Log-Likelihood
#'
#' Computes the observed-data log-likelihood for the multiway mixture model,
#' given the observed data, mean estimates, variance estimates, mixing proportions,
#' and partially observed labels.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param M Numeric array of mean estimates with dimensions
#'   \code{(p x C x V)}, where \code{C} is the number of conditions
#'   and \code{V} is the number of infection status
#' @param sigma2 Numeric array of variance estimates with dimensions
#'   \code{(p x C x V)}.
#' @param probs Numeric matrix of mixing proportions with dimensions
#'   \code{(C x V)}.
#' @param c_obs Integer vector of length \code{n}, giving the observed
#'   condition labels for each sample (use \code{NA} if unknown).
#' @param v_obs Integer vector of length \code{n}, giving the observed
#'   variant labels for each sample (use \code{NA} if unknown).
#'
#' @return A numeric scalar: the observed-data log-likelihood value.
#'
#' @importFrom stats dnorm
#' @importFrom matrixStats logSumExp
#' @keywords internal
observed.data.loglik <- function(Y, M, sigma2, probs, c_obs, v_obs){

  n <- dim(Y)[1]; p <- dim(Y)[2]
  c_dim <- dim(M)[2]; v_dim <- dim(M)[3]
  cond_probs_row <- probs/rowSums(probs)
  cond_probs_col <- sweep(probs, 2, colSums(probs), FUN = "/")
  tot.contrib.vec <- rep(0, n)
  
  tot.contrib <- NULL
  both_known_contrib <- 0
  v_unknown_contrib <- 0
  c_unknown_contrib <- 0
  both_unknown_contrib <- 0
  
  # V observed and C observed
  both_known <- which(!is.na(v_obs) & !is.na(c_obs)) 
  if (length(both_known) > 0) {
    both_known_contrib <- sapply(both_known, function(index) {
        p0 <- probs[c_obs[index], v_obs[index]]
        if (is.na(p0) | p0 == 0) p0 <- .001
        sum(dnorm(Y[index, ], mean = M[, c_obs[index], v_obs[index]],
                sd = sqrt(sigma2[, c_obs[index], v_obs[index]]),
                log = TRUE)) + log(p0)
    })
  }
  
  # C unknown
  c_unknown <- which(!is.na(v_obs) & is.na(c_obs))
  if (length(c_unknown) > 0) {
    combos <- expand.grid(i = c_unknown, cc = 1:c_dim)
    combos$c_unknown_contrib <- mapply(function(i, cc) {
      p1 <- cond_probs_col[cc, v_obs[i]]
      if (is.na(p1) | p1 == 0) p1 <- .001
      logSumExp(dnorm(Y[i, ], mean = M[, cc, v_obs[i]],
                      sd = sqrt(sigma2[, cc, v_obs[i]]), log = TRUE)) +
        log(p1)}, combos$i, combos$cc)
    c_unknown_contrib <- tapply(combos$c_unknown_contrib, combos$i, logSumExp)
  }
  
  # V unknown
  v_unknown <- which(is.na(v_obs) & !is.na(c_obs))
  if (length(v_unknown) > 0) {
    combos <- expand.grid(i = v_unknown, v = 1:v_dim)
    combos$v_unknown_contrib <- mapply(function(i, v) {
        p2 <- cond_probs_row[c_obs[i], v]
        if(is.na(p2) | p2 == 0) p2 <- .001
      sum(dnorm(Y[i, ], mean = M[, c_obs[i], v],
                sd = sqrt(sigma2[, c_obs[i], v]), log = TRUE)) + log(p2)
        }, combos$i, combos$v)
    v_unknown_contrib <- tapply(combos$v_unknown_contrib, combos$i, logSumExp)
  }
  
  # Both unknown
  both_unknown <- which(is.na(v_obs) & is.na(c_obs))
  if (length(both_unknown) > 0) {
    combos <- expand.grid(i = both_unknown, v = 1:v_dim, cc=1:c_dim)
    combos$both_unknown_contrib <- mapply(function(i, v, cc) {
        p3 <- probs[cc, v]
        if(is.na(p3) | p3 == 0) p3 <- .001
      logSumExp(dnorm(Y[i, ], mean = M[, cc, v],
                      sd = sqrt(sigma2[, cc, v]), log = TRUE)) +
        log(p3)}, combos$i, combos$v, combos$cc)
    both_unknown_contrib <- tapply(combos$both_unknown_contrib, combos$i, logSumExp)
  }
  
  # Sum together
  tot.contrib <- sum(c(both_unknown_contrib, v_unknown_contrib, 
                       c_unknown_contrib, both_known_contrib))
  
  return(tot.contrib)
}
