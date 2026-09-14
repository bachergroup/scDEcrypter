#' Raw Seurat Object with Single-Cell RNA-Seq Data
#'
#' A Seurat object containing raw single-cell RNA-seq count data. This is a subset of 
#' 2,000 cells (of MOCK and DPI3) from Ravindra et al., 2021.
#' 
#' @format A Seurat object with gene expression counts stored in the RNA assay.
#'   The object contains raw (non-normalized) count data for cells and genes.
#'
#' @details
#' This is a raw Seurat object that can be preprocessed using the 
#' \code{\link{preprocess_scDEcrypter}} function for differential expression analysis.
#' The object includes all standard Seurat metadata and assays.
#'
#' @examples
#' \dontrun{
#' data(seu_sub)
#' # View the object
#' seu_sub
#' }
"seu_sub"
