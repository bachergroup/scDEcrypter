
#' Classification Accuracy
#'
#' Computes the proportion of correctly predicted labels compared
#' to the true labels in a classification task.
#'
#' @param indat A data frame containing at least two columns:
#'   \itemize{
#'     \item \code{Prediction} — predicted class labels.
#'     \item \code{True_infection} — true class labels.
#'   }
#'
#' @return A numeric scalar between 0 and 1, representing the
#'   classification accuracy.
#'
#' @export
acc1 <- function(indat) {
  print(nrow(indat))
  sum(indat$Prediction == indat$True_infection) / nrow(indat)
}

