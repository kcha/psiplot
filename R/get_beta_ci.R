#' Calculate confidence intervals for error bars
#'
#' Helper function to filter and return confidence intervals based on beta
#' distribution from Q scores. For internal use. Beta sampling functions
#' provided by Tim Sterne-Weiler.
#'
#' @param q a data frame of PSI and corresponding quality values
#' @return Confidence intevals of PSI values
#' @author Tim Sterne-Weiler, Kevin Ha
get_beta_ci <- function(q) {
  #  This function takes a qual and returns c(post_alpha, post_beta)
  #  Increments by prior alpha and prior distribution beta, uniform by default
  parseQual <- function(qual, prior_alpha=1, prior_beta=1) {
    if(is.na(qual) || !grepl("@", qual)) { return(c(prior_alpha, prior_beta)) }  ## for INT NA Columns
    res <- unlist(strsplit(unlist(strsplit(as.character(qual), "@"))[2], ","))
    if(is.na(res[1]) || is.na(res[2])) { return(c(prior_alpha, prior_beta)) }
    if(is.null(res[1]) || is.null(res[2])) { return(c(prior_alpha, prior_beta)) }
    if(res[1] == "NA" || res[2] == "NA") { return(c(prior_alpha, prior_beta)) }
    res <- as.numeric(res)
    if(is.nan(res[1]) || is.nan(res[2])) { return(c(prior_alpha, prior_beta)) }
    if(is.infinite(res[1]) || is.infinite(res[2])) { return(c(prior_alpha, prior_beta)) }
    res[1] <- res[1] + prior_alpha
    res[2] <- res[2] + prior_beta
    res
  }

  betaCI <- function(betaDist, percentile = c(0.05, 0.95)) {
    quantile(betaDist, p=percentile, na.rm = T)
  }

  # Extention of betaCI function that includes the sampling step
  betaCISample <- function(alpha, beta, n = 5000) {
    if (is.na(alpha) || is.na(beta)) {
      sample <- NA
    } else {
      set.seed(79)
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

