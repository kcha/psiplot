#' Format table to keep only PSIs and convert exon metadata as rownames.
#'
#' Calls \code{convert_psi()}.
#'
#' @param m Original PSI plus quality score data frame
#' @return Formatted data frame keeping only PSI values
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
