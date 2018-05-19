#' Convert low/bad quality PSI values to NA
#'
#' Helper function to filter and return PSI values.
#' PSIs are converted to NA if first coverage code is 'N'
#' e.g. PSI=100, Coverage=N,N,N,OK,S ---> PSI=NA
#' For internal use. Called by \code{\link{format_table}}.
#'
#' @param t Original PSI plus quality scores data frame WITHOUT the exon
#' metadata columns
#' @return Data frame with the same dimensions as \emph{t} and low/bad quality P
#' SI values converted to \code{NA}
#' @examples
#' \dontrun{
#' psiplot:::convert_psi(psi[,7:ncol(psi)])
#' }
convert_psi <- function(t) {
  stopifnot(ncol(t) %% 2 == 0)
  psi <- t

  for (i in seq(1, ncol(psi), 2)) {
    cov <- strsplit(as.character(psi[,i+1]), split = ",")
    cov <- sapply(cov, "[", 1)

    na <- which(cov == "N")
    if (length(na) > 0)
      psi[na, i] <- NA
  }
  return(psi)
}
