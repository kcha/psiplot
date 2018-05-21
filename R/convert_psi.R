#' Convert low/bad quality PSI values to NA
#'
#' Helper function to filter and return PSI values.
#' PSIs are converted to NA if first coverage code is 'N'
#' e.g. PSI=100, Coverage=N,N,N,OK,S ---> PSI=NA
#' For internal use. Called by \code{\link{format_table}}.
#'
#' @param t Original PSI plus quality scores data frame WITHOUT the exon
#' metadata columns
#' @param qual String indicating the minimun \emph{vast-tools} quality score
#' for the PSI to be accepted. Defaults to \code{'VLOW'}. See the
#' \href{https://github.com/vastgroup/vast-tools/blob/master/README.md}{vast-tools
#' documentation} for details.
#' @return Data frame with the same dimensions as \emph{t} and low/bad quality P
#' SI values converted to \code{NA}
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @examples
#' \dontrun{
#' psiplot:::convert_psi(psi[,7:ncol(psi)])
#' }
convert_psi <- function(t,qual = c("VLOW","N","LOW","OK","SOK")) {
  qual <- match.arg(qual)
  stopifnot(ncol(t) %% 2 == 0)
  psi <- t

  qual_vector <- c("N","VLOW","LOW","OK","SOK")
  min_qual <- which(qual_vector==qual)

  for (i in seq(1, ncol(psi), 2)) {
    cov <- apply(psi[,i+1],
                 1,
                 function(x) unlist(strsplit(x,split=",",fixed = T))[[1]])
    cov <- match(cov,qual_vector)
    na <- which(cov < min_qual)
    if (length(na) > 0)
      psi[na, i] <- NA
  }
  return(psi)
}
