#' Format input event data
#' 
#' Convert low/bad quality PSI values and convert event metadata as rownames.
#' Prepares input event data for plotting. Calls \code{\link{convert_psi}}.
#'
#' @param m A data frame of event data as outputted by
#' \code{vast-tools combine}. e.g. each row is an event containing exon 
#' metadata, PSI and quality scores values.
#' @return Data frame with PSI and quality scores. Rownames set as a
#' concatenation of exon metadata delimited by |.
#' @export
#' @examples
#' # For example input, see:
#' mm.psi
#' format_table(mm.psi)
format_table <- function(m) {
  if (!grepl("^GENE", colnames(m)[1])) {
    stop("Invalid column names. Does your input file contain the correct header?")
  }
  # Format table to keep only PSIs and convert exon metadata as rownames
  id <- paste(m$COMPLEX, m$GENE, m$COORD, m$LENGTH, sep="|")

  # Extract PSIs
  r <- convert_psi(m[,7:ncol(m)])
  rownames(r) <- id
  return(r)
}
