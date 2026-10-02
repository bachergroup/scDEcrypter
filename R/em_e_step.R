#' E-step of Multiway Mixture Model
#'
#' @param Y Expression matrix (cell x gene).
#' @param c_obs Observed condition labels (vector).
#' @param v_obs Observed infection status labels (vector).
#' @param M Mean array (gene x condition x infection status).
#' @param probs Mixing proportions (condition x infection status matrix).
#' @param sigma2 Variance array (gene x condition x infection status).
#' @return Weight tensor (cell x condition x infection status).
#' @keywords internal
E_step <- function(Y, c_obs, v_obs, M, probs, sigma2) {
    Y <- as.matrix(Y)
    if (length(dim(Y)) != 2) {
        stop("Y must be a 2D matrix with cells in rows and genes in columns.", call. = FALSE)
    }

    c_obs <- as.integer(c_obs)
    v_obs <- as.integer(v_obs)

    W <- E_step_rcpp(Y, c_obs, v_obs, M, probs, sigma2)

    y_rownames <- rownames(Y)
    if (!is.null(y_rownames) && length(y_rownames) == nrow(Y)) {
        dimnames(W)[[1]] <- y_rownames
    }
    if (!is.null(dimnames(M)) && length(dimnames(M)) >= 3) {
        dimnames(W)[2:3] <- dimnames(M)[2:3]
    }

    W
}