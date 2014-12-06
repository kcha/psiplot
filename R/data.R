#' Sample genes with simulated PSI data
#'
#' Contains simulated PSI and quality score data as produced by \code{vast-tools
#' combine} of a few randomly selected mouse genes.
#' 
#' @docType data
#' @name mm.psi
#' @usage mm.psi
#' @format A 5 x 14 data frame
#' @keywords datasets
NULL

#' Sample psiplot configuration settings for dataset \code{mm.psi}
#'
#' Example of how a psiplot configuration file should be formatted. This can be
#' passed to \code{psiplot::plot_event} using the \code{config} argument, which
#' accepts either the file path of the config file or an \code{n x 4} data frame where
#' \code{n} is the number of samples.
#' @docType data
#' @name mm.psi.config
#' @usage mm.psi.config
#' @format A 4 x 4 data frame
#' @keywords datasets
NULL