
#' Gene- and Cell-Type-Specific Observed Data Log-Likelihood
#'
#' Computes the observed data log-likelihood for each gene within each
#' cell type, given the estimated mixture parameters and observed
#' conditions/contrasts.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param M A numeric array of dimension \eqn{(p \times C \times V)},
#'   representing the estimated mean expression values for each
#'   gene \eqn{p}, cell type \eqn{C}, and condition/contrast \eqn{V}.
#' @param sigma2 A numeric array of the same dimension as \code{M},
#'   representing the estimated variances.
#' @param probs A numeric matrix of dimension \eqn{(C \times V)},
#'   containing the mixture probabilities for each cell type and
#'   condition.
#' @param Condition An integer or factor vector of length \eqn{n},
#'   specifying the observed cell type labels for each sample.
#'   Missing values (\code{NA}) indicate unobserved cell types.
#' @param Contrast An integer or factor vector of length \eqn{n},
#'   specifying the observed condition/contrast labels for each sample.
#'   Missing values (\code{NA}) indicate unobserved conditions.
#'
#' @return A numeric matrix of dimension \eqn{(p \times C)}, where each
#'   entry \code{LL_gene_celltype[gene, celltype]} gives the
#'   log-likelihood contribution of that gene under the given
#'   cell type.
#'
#' @export
observed.data.loglik.gene_celltype <- function (Y, M, sigma2, probs, Condition, Contrast)
{
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  C.dim <- length(unique(Condition))
  V.dim <- length(unique(Contrast))
  LL_gene_celltype <- matrix(0, nrow = p, ncol = C.dim)
  cond.probs.row <- probs/rowSums(probs)
  cond.probs.col <- sweep(probs, 2, colSums(probs), FUN = "/")

  for (c in 1:C.dim) {

    both_known <- which(!is.na(Contrast) & !is.na(Condition) & Condition == c)
    tot.contrib_1 <- array(0, dim=c(p))
    log_likelihood_kc <- NULL
    for(kk in 1:p){
      tot.contrib_1[kk] <- sum(dnorm(Y[both_known, kk],
                                     mean = M[kk, c, Contrast[both_known]],
                                     sd = sqrt(sigma2[kk, c, Contrast[both_known]]),
                                     log = TRUE))
    }
    log_likelihood_kc <- tot.contrib_1

    if (any(!is.na(Contrast) & is.na(Condition))) {
      v_known <- which(!is.na(Contrast) & is.na(Condition))
      tot.contrib_2 <- array(0, dim=c(p))
      for(kk in 1:p){
        tot.contrib_2[kk] <- sum(dnorm(Y[v_known, kk],
                                       mean = M[kk, c, Contrast[v_known]],
                                       sd = sqrt(sigma2[kk, c, Contrast[v_known]]),
                                       log = TRUE) + log(cond.probs.col[c, Contrast[v_known]]))
      }
      log_likelihood_kc <- tot.contrib_1 + tot.contrib_2
    }

    LL_gene_celltype[, c] <- log_likelihood_kc
  }

  return(LL_gene_celltype)
}
