
#' Split Data into Training and Test Sets
#'
#' Returns the indices to split the observed dataset into training 
#' ("generation") and test subsets, stratified by labels. 
#' Cells with missing variant labels are also split proportionally.
#'
#' @param c_obs Integer vector of length \code{n}, giving the pre-labeled
#'   partitioning variable status for each cell (may contain \code{NA}).
#' @param v_obs Integer vector of length \code{n}, giving the pre-labeled
#'   infection status labels for each cell (may contain \code{NA}).
#' @param train_frac Numeric in \code{(0,1)}. Proportion of cells
#'   within each variant (and missing group) assigned to the training
#'   set. Defaults to \code{0.5}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{generation_idx}}{Training indices.}
#'   \item{\code{test_idx}}{Testing indices.}
#' }
#'
#' @importFrom stats na.omit
splitDataIDX <- function(c_obs, v_obs, train_frac = 0.5) {
    
    v_levels <- sort(unique(na.omit(v_obs)))
    
    generation_idx <- c()
    test_idx <- c()
    for (v in v_levels) {
        idx_v <- names(which(v_obs == v))
        n_v <- length(idx_v)
        n_train <- round(train_frac * n_v)
        train_idx_v <- sample(idx_v, n_train)
        test_idx_v <- setdiff(idx_v, train_idx_v)
        
        generation_idx <- c(generation_idx, train_idx_v)
        test_idx <- c(test_idx, test_idx_v)
    }
    
    idx_na <- names(which(is.na(v_obs)))
    if (length(idx_na) > 0) {
        n_na <- length(idx_na)
        n_train_na <- round(train_frac * n_na)
        train_idx_na <- sample(idx_na, n_train_na)
        test_idx_na <- setdiff(idx_na, train_idx_na)
        
        generation_idx <- c(generation_idx, train_idx_na)
        test_idx <- c(test_idx, test_idx_na)
    }
    
    return(list(generation_idx = generation_idx, test_idx=test_idx))
}



#' Preprocess data for scDEcrypter
#' 
#' Data will first be split into training and testing, followed
#' by independent variance stabilizing normalization. 
#' 
#' @param Data A Seurat object containing an RNA assay with raw counts data.
#' @param c_obs_var Character name of meta data column giving the pre-labeled
#'   partitioning variable status for each cell (may contain \code{NA}).
#' @param v_obs_var Character name of meta data column giving the pre-labeled
#'   infection status labels for each cell (may contain \code{NA}).
#' @param train_frac Numeric in \code{(0,1)}. Proportion of cells
#'   within each variant (and missing group) assigned to the training
#'   set. Defaults to \code{0.5}.
#' @param vs_method Transformation method passed to transformGamPoi. Default is \code{shifted_log}.
#' @importFrom transformGamPoi transformGamPoi
#' @import Seurat
#' @export
preprocess_scDEcrypter <- function(Data, 
                                   c_obs_var = "C.preLabel", 
                                   v_obs_var = "V.preLabel", 
                                   train_frac = .5, 
                                   vs_method = c("shifted_log", "none")) {
    
    vs_method <- match.arg(vs_method)
    
    if (!"RNA" %in% names(Data@assays)) {
        stop("Seurat object must contain an RNA assay.")
    }
    if (!all(c(c_obs_var, v_obs_var) %in% colnames(Data@meta.data))) {
      stop("One of the pre-labeled columns (named in c_obs_var v_obs_var) is not found
             in the Seurat metadata. Please check spelling and try again.")
    }
    message("Splitting data into generation and test sets...")
    index_splits <- splitDataIDX(Data[[c_obs_var, drop=TRUE]], 
                                 Data[[v_obs_var, drop=TRUE]], 
                                 train_frac = train_frac)
    
    Data$Index_Split <- ""
    Data$Index_Split[index_splits$generation_idx] <- "Generation"
    Data$Index_Split[index_splits$test_idx] <- "Test"
    
    Data[["RNA"]] <- split(Data[["RNA"]], f = Data$Index_Split)

    message("Applying VST to generation set...")    
    Y_gen <- transformGamPoi::transformGamPoi(as.matrix(Data[["RNA"]]$counts.Generation),
                                                   transformation = vs_method,
                                                   size_factors = "poscounts"
    )
    message("Applying VST to test set...")        
    Y_test <- transformGamPoi::transformGamPoi(as.matrix(Data[["RNA"]]$counts.Test),
                                                    transformation = vs_method,
                                                    size_factors = "poscounts"
    )
    
    Data[["RNA"]]$data.Generation <- Y_gen
    Data[["RNA"]]$data.Test <- Y_test
    
    return(Data)
}


#' Identify Top Variable Genes over Cell Partitions
#'
#' Selects the most variable genes within each group of cells defined by a partitioning column in a Seurat object.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param partition_colname Character. Name of the cell-level metadata column used to partition cells into groups (e.g., cell type, cluster).
#' @param n_genes_per_group Integer. Number of top variable genes to rank per group (default: 500).
#'
#' @return A character vector of unique gene names with highest variance across the specified groups.
#'
#' @details
#' For each group defined by \code{partition_colname}, the function calculates the variance of each gene across cells in that group,
#' selects the top \code{n_genes_per_group} most variable genes, and then returns the unique set of the top \code{n_top} genes from each group.
#'
#' @importFrom stats setNames var
#' @export

get_top_genes <- function(seurat_obj, partition_colname="C.preLabel", 
													n_genes_per_group = 50) {

  vargroup <- seurat_obj[[partition_colname]]
	vargroup <- setNames(vargroup[[1]], rownames(vargroup))
	
  top_genesA <- sapply(unique(vargroup), function(x) {
							tmpCells <- names(which(vargroup == x))
							all_var <- apply(seurat_obj[["RNA"]]$data.Generation[, tmpCells], 1, var)
							all_var_ord <- all_var[order(all_var, decreasing = TRUE)[1:n_genes_per_group]]
    return(names(all_var_ord))
  })
  top_genes <- unique(as.vector(top_genesA))
  return(top_genes)
}


