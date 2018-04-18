#' Get samples
#'
#' Get samples from PSI input table
#'
#' @param df A data frame of PSI input values.
#' @param configdf A 4-column data frame.
#' @param value If \code{TRUE}, then returns the sample names. Otherwise, return
#' the indices.
#' @return If \code{value = TRUE}, A character vector consisting of sample names.
#' Otherwise, a numeric vector of indices.
#' @examples
#' get_psi_samples(psi)
get_psi_samples <- function(df, configdf = NULL, value = TRUE) {
  ix <- seq(7, ncol(df), 2)
  samples <- colnames(df)[ix]
  if (!is.null(configdf)) {
    if (!is.data.frame(configdf)) {
      stop("configdf is not a data frame!")
    }

    samples <- configdf$SampleName[which(configdf$SampleName %in% samples)]
  }
  if (value) {
    return(samples)
  }
  return(match(samples, colnames(df)))
}
