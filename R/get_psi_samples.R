#' Get samples from PSI input table
#'
#' @details
#'
#' When \code{config} is provided, then the samples from the \emph{SampleName}
#' column are used to find the indices in the PSI table.
#'
#' @param df A data frame of PSI input values.
#' @param config Config data frame
#' @param value If \code{TRUE}, then returns the sample names. Otherwise, return
#' the indices.
#' @return If \code{value = TRUE}, A character vector consisting of sample names.
#' Otherwise, a numeric vector of indices.
#' @export
#' @examples
#' \dontrun{
#' get_psi_samples(psi)
#' }
get_psi_samples <- function(df, config = NULL, value = TRUE) {
  sample_cols <- colnames(df[, 7:ncol(df)])
  samples <- grep("(\\.|-)Q$", sample_cols, invert = TRUE, value = TRUE)
  # samples <- colnames(df)[ix]
  if (!is.null(config)) {
    if (!is.data.frame(config)) {
      stop("configdf is not a data frame!")
    }

    samples <- config$SampleName[which(config$SampleName %in% samples)]
  }
  if (value) {
    return(samples)
  }
  return(match(samples, colnames(df)))
}
