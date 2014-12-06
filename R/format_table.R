#' Format input event data
#' 
#' Format table to keep only PSIs and convert exon metadata as rownames.
#' Prepares input event data for plotting. Calls \code{convert_psi()}.
#'
#' @param m A data frame of PSI and quality scores as outputted by
#' \code{vast-tools combine}
#' @return Formatted data frame keeping only PSI values
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
