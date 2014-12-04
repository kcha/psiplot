#' Convert low/bad quality PSI values to NA
#'
#' Helper function to filter and return PSI values
#' PSIs are converted to NA if first coverage code is 'N'
#' e.g. PSI=100, Coverage=N,N,N,OK,S ---> PSI=NA
#'
#' @param t Original PSI plus quality scores data frame WITHOUT the first 7 columns
#' @return Data frame with converted NA PSI values
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
