

#' Multiway Mixture Model with Penalty
#'
#' Fits a multiway mixture model using the EM algorithm with
#' an \eqn{\ell_2}-type penalty on mean parameters to encourage
#' similarity across variants.
#'
#' @param Data Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param c_obs Integer vector of length \code{n}, giving the pre-labeled
#'   partitioning variable status for each cell (may contain \code{NA}).
#' @param v_obs Integer vector of length \code{n}, giving the pre-labeled
#'   infection status labels for each cell (may contain \code{NA}).
#' @param max.iter Integer. Maximum number of EM iterations for each
#'   value of \code{lambda}.
#' @param tol Numeric. Convergence tolerance for relative change in
#'   mean estimates.
#' @param c_star Integer. Number of condition clusters to fit.
#' @param v_star Integer. Number of variant clusters to fit.
#' @param lambda.vec Numeric vector of penalty parameters. A separate
#'   model is fitted for each \code{lambda}.
#' @param infectionLabels Character vector of length equal to the number of viral types,
#'   specifying biological names for each infection/variant status. 
#'   Example: \code{c("Uninfected", "Infected", "Bystander", etc)}
#' @param partitionLabels Character vector of length equal to the number of cell types,
#'   specifying biological names for each cell type cluster.
#'   Example: \code{c("T_cells", "B_cells", "Monocytes", etc)}
#' @return A list containing:
#' \describe{
#'   \item{\code{M_list}}{List of estimated means (\code{p x C x V})
#'   for each \code{lambda}.}
#'   \item{\code{sigma2_list}}{List of estimated variances
#'   (\code{p x C x V}) for each \code{lambda}.}
#'   \item{\code{probs_list}}{List of estimated mixing proportions
#'   (\code{C x V}) for each \code{lambda}.}
#'   \item{\code{weights_list}}{List of posterior responsibilities
#'   (\code{n x C x V}) for each \code{lambda}.}
#' }
#'
#' @export
fit_scDEcrypter <- function(Data, c_obs=NULL, v_obs=NULL, infectionLabels=NULL, partitionLabels=NULL,
                      max.iter=250, tol=1e-8, c_star, v_star, lambda.vec){


  if (!"RNA" %in% names(Data@assays)) {
      stop("Seurat object must contain an RNA assay.")
  }
 
  if(is.null(c_obs)) c_obs <- Data$C.preLabel
  if(is.null(v_obs)) v_obs <- Data$V.preLabel

  # Extract and rotate data
  Y <- t(Data[["RNA"]]$data.Generation)

  M_lambda_list <- list()
  sigma2_lambda_list <- list()
  probs_lambda_list <- list()
  weights_lambda_list <- list()
  message("Initializing...")
  tmp <- initializer_scDEcrypter(Y, c_obs, v_obs, max.iter, tol, c_star, v_star)
  
  
  for(lambda in lambda.vec) {
    
    M <- tmp$M
    sigma2 <- tmp$sigma2
    probs <- tmp$probs
    
    message("EM Iterating...")
    for(mm in seq_len(max.iter)){
      W1 <- E_step(Y, c_obs, v_obs, M, probs, sigma2)
      probs.new <- M_step_probs(Y, W1)
      M.new <- update_mu(Y, M, sigma2, W1, lambda = lambda * nrow(Y))
      sigma2.new <- M_step_variance(Y, W1, M.new)
      
      cat(sum((M - M.new)^2)/sum(M^2), "\n")
      if(sum((M - M.new)^2)/sum(M^2) < tol){
        break
      }
      M <- M.new
      sigma2 <- sigma2.new
      probs <- probs.new
    }
    message("Preparing outputs...")
    out.weights <- E_step(Y, c_obs, v_obs, M, probs, sigma2)
    
    dimnames(M.new) <- list(dimnames(M.new)[[1]], 
                               paste0("partitionStatus_", seq_len(c_star)),
                               paste0("infectionStatus_", seq_len(v_star)))
    dimnames(sigma2.new) <- dimnames(M.new)
    dimnames(out.weights) <- list(dimnames(out.weights)[[1]],
                                 paste0("partitionStatus_", seq_len(c_star)),
                                 paste0("infectionStatus_", seq_len(v_star)))
    dimnames(probs.new) <- list(paste0("partitionStatus_", seq_len(c_star)),
                                 paste0("infectionStatus_", seq_len(v_star)))
                                 
    M_lambda_list[[as.character(lambda)]] <- M.new
    sigma2_lambda_list[[as.character(lambda)]] <- sigma2.new
    probs_lambda_list[[as.character(lambda)]] <- probs.new
    weights_lambda_list[[as.character(lambda)]] <- out.weights
  }
  
  if (length(lambda.vec) == 1) {
    out.list <- list(
      "M_generation" = M_lambda_list[[1]],
      "sigma2_generation" = sigma2_lambda_list[[1]], 
      "probs_generation" = probs_lambda_list[[1]],
      "weights_generation" = weights_lambda_list[[1]]
    )
    if(!is.null(infectionLabels) & !is.null(partitionLabels)) {
     out.list <- switchNames(out.list, infectionLabels, partitionLabels)
    }
    return(out.list)
  } else {                             
    return(list(
      "M_generation" = M_lambda_list,
      "sigma2_generation" = sigma2_lambda_list, 
      "probs_generation" = probs_lambda_list,
      "weights_generation" = weights_lambda_list
    ))
  }

}