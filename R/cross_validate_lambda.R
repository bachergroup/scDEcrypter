
#' Cross-validation for Lambda Selection
#'
#' Performs K-fold cross-validation to select the best penalty parameter
#' \eqn{\lambda} for the \code{MultiwayMixture} model. Each fold splits the
#' data into training and validation sets, fits the model on training data
#' across candidate \eqn{\lambda} values, and evaluates log-likelihood on
#' validation data. The best \eqn{\lambda} is chosen by maximizing the
#' average validation log-likelihood across folds.
#'
#' @param Y_generation Numeric matrix of dimension \code{n x p}, training data.
#'   Each row is a sample (cell), each column is a feature (e.g. gene).
#' @param C.obs_generation Integer vector of length \code{n}, observed condition
#'   labels (may contain \code{NA}).
#' @param V.obs_generation Integer vector of length \code{n}, observed variant
#'   labels (may contain \code{NA}).
#' @param lambda.vec Numeric vector. Candidate \eqn{\lambda} values to be
#'   cross-validated.
#' @param max.iter Integer. Maximum number of EM iterations for
#'   \code{MultiwayMixture}.
#' @param tol Numeric. Convergence tolerance for EM updates.
#' @param C.star Integer. Prespecified number of latent conditions.
#' @param V.star Integer. Prespecified number of latent variants.
#' @param seed Integer. Random seed for reproducibility.
#' @param NCORES Integer. Number of cores for parallel computation.
#' @param K Integer. Number of cross-validation folds. Defaults to 5.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{best_lambda}}{The \eqn{\lambda} value with highest average
#'   validation log-likelihood.}
#'   \item{\code{lambda_log_likelihoods}}{Matrix of log-likelihoods, rows =
#'   \code{lambda.vec}, columns = folds.}
#'   \item{\code{avg_log_likelihoods}}{Vector of mean log-likelihoods across folds.}
#' }
#' @export

cross_validate_lambda <- function(Y_generation, C.obs_generation, V.obs_generation,
                                  lambda.vec, max.iter, tol, C.star, V.star, seed, NCORES, K = 5) {
  
	set.seed(seed)
	print("Starting CrossValidation")

	V.levels <- sort(unique(na.omit(V.obs_generation)))
	folds_list <- vector("list", K)

	for (v in V.levels) {
	  idx_v <- which(V.obs_generation == v)
	  fold_indices_v <- createFolds(1:length(idx_v), k = K, list = TRUE, returnTrain = FALSE)
	  for (k in 1:K) {
	    folds_list[[k]] <- c(folds_list[[k]], idx_v[fold_indices_v[[k]]])
	  }
	}

	idx_na <- which(is.na(V.obs_generation))
	if (length(idx_na) > 0) {
	  fold_indices_na <- createFolds(1:length(idx_na), k = K, list = TRUE, returnTrain = FALSE)
	  for (k in 1:K) {
	    folds_list[[k]] <- c(folds_list[[k]], idx_na[fold_indices_na[[k]]])
	  }
	}

	lambda_log_likelihoods <- mclapply(1:K, runFOLDS, folds_list, nrow(Y_generation),
	                                   Y_generation, C.obs_generation, V.obs_generation,
	                                   max.iter, tol, C.star, V.star, seed, lambda.vec,
	                                   mc.cores = NCORES)

	lambda_log_likelihoods <- do.call(cbind, lambda_log_likelihoods)
	colnames(lambda_log_likelihoods) <- paste0("Fold_", 1:K)
	rownames(lambda_log_likelihoods) <- lambda.vec

	avg_log_likelihoods <- rowMeans(lambda_log_likelihoods)
	best_lambda <- lambda.vec[which.max(avg_log_likelihoods)]
	cat("Best lambda selected: ", best_lambda, "\n")

	return(list(
	  "best_lambda" = best_lambda,
	  "lambda_log_likelihoods" = lambda_log_likelihoods,
	  "avg_log_likelihoods" = avg_log_likelihoods
	))
}

#' @keywords internal
runFOLDS <- function(x, fold_indices, n, Y_generation, C.obs_generation,
                     V.obs_generation, max.iter, tol, C.star, V.star,
                     seed, lambda.vec) {
											 
	set.seed(seed)
	validation_idx <- fold_indices[[x]]
	training_idx <- setdiff(1:n, validation_idx)

	Y_train <- Y_generation[training_idx, ]
	C.obs_train <- C.obs_generation[training_idx]
	V.obs_train <- V.obs_generation[training_idx]

	Y_val <- Y_generation[validation_idx, ]
	C.obs_val <- C.obs_generation[validation_idx]
	V.obs_val <- V.obs_generation[validation_idx]

	print(paste0("Running Folds ", x))
	# Perform lambda selection on training set
	results <- MultiwayMixture(Y = Y_train, C.obs = C.obs_train, V.obs = V.obs_train,
	                          max.iter = max.iter, tol = tol, C.star = C.star, V.star = V.star,
	                          seed = seed, lambda.vec = length(training_idx)*lambda.vec)
	print(paste0("Done Folds ", x))
	# Evaluate each lambda on the validation set
	log_likelihood <- NULL
	for (lambda_idx in 1:length(lambda.vec)) {
	 lambda <- lambda.vec[lambda_idx]
	 M_lambda <- results$M_list[[lambda_idx]]
	 sigma2_lambda <- results$sigma2_list[[lambda_idx]]
	 probs_lambda <- results$probs_list[[lambda_idx]]

	 # Calculate log-likelihood for this lambda on the validation set
	 log_likelihood <- c(log_likelihood, observed.data.loglik(Y_val, M_lambda, sigma2_lambda, probs_lambda, C.obs_val, V.obs_val))
	}
	print(paste0("Done likelihood calculation ", x))
	return(log_likelihood)
}
