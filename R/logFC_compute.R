
#' Compute Log-Fold Changes Between Conditions
#'
#' Computes log-fold change (logFC) contrasts between all pairs of
#' conditions/contrasts for each gene and cell type, given estimated
#' mean expression values.
#'
#' @param M A numeric array of dimension \eqn{(p \times C \times V)},
#'   where \eqn{p} is the number of genes, \eqn{C} is the number of
#'   cell types, and \eqn{V} is the number of conditions. Entries
#'   represent estimated mean expression levels.
#' @param Condition Optional character vector of length \eqn{C},
#'   specifying names of cell types. If \code{NULL}, generic labels
#'   \code{"celltype1"}, \code{"celltype2"}, … are used in the output.
#' @param V.levels Optional character vector of length \eqn{V},
#'   specifying names of conditions. If \code{NULL}, defaults to
#'   \code{"V1"}, \code{"V2"}, ….
#'
#' @return A named list of data frames, one for each condition contrast
#'   (e.g., \code{"V1_vs_V2"}). Each data frame has \eqn{p} rows
#'   (genes) and \eqn{C} columns (cell types), where each entry is the
#'   log-fold change between the two conditions:
#'   \deqn{ \text{logFC}_{j,c} = M_{j,c,v1} - M_{j,c,v2} }
#'
#' @export
logFC_compute <- function(M, Condition = NULL, V.levels = NULL) {
  num_genes <- dim(M)[1]
  num_celltypes <- dim(M)[2]
  V_dim <- dim(M)[3]

  if (is.null(V.levels)) {
    V.levels <- paste0("V", 1:V_dim)
  }

  combs <- combn(V_dim, 2)
  logFC_list <- list()

  for (i in 1:ncol(combs)) {
    v1 <- combs[1, i]
    v2 <- combs[2, i]

    df <- matrix(NA, nrow = num_genes, ncol = num_celltypes)
    for (c in 1:num_celltypes) {
      df[, c] <- M[, c, v1] - M[, c, v2]
    }

    colnames(df) <- if (!is.null(Condition)) {
      paste0("logFC_", Condition, "_", V.levels[v1], "_vs_", V.levels[v2])
    } else {
      paste0("logFC_celltype", 1:num_celltypes, "_", V.levels[v1], "_vs_", V.levels[v2])
    }

    logFC_list[[paste0(V.levels[v1], "_vs_", V.levels[v2])]] <- as.data.frame(df)
  }

  return(logFC_list)
}


