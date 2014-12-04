#' Calculate confidence intervals for error bars
#'
#' Helper function to filter and return confidence intervals based on beta
#' distribution from Q scores. For internal use.
#'
#' @param q original PSI plus quality scores WITHOUT the first 7 columns
#' @return Confidence intevals of PSI values
get_beta_ci <- function(q) {
  betaCI <- function(betaDist, percentile = c(0.05, 0.95)) {
    quantile(betaDist, p=percentile, na.rm = T)
  }

  # Extention of betaCI function that includes the sampling step
  betaCISample <- function(alpha, beta, n = 5000) {
    if (is.na(alpha) || is.na(beta)) {
      sample <- NA
    } else {
      sample <- rbeta(n, alpha, beta)
    }
    return(betaCI(sample))
  }

  parameters <- sapply(q, function(x) parseQual(as.character(x)))
  ci <- lapply(1:ncol(parameters), function(j) betaCISample(parameters[1,j],
                                                            parameters[2,j]))
  ci <- do.call("rbind", ci) * 100
  return(ci)
}

