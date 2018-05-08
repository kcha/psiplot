#' Format input event data
#'
#' Convert low/bad quality PSI values and convert event metadata as rownames.
#' Prepares input event data for plotting. Calls \code{\link{convert_psi}}.
#'
#' @details
#' This function is also used for formatting cRPKM input data by setting
#' \code{expr = TRUE}.
#'
#' @param x A data frame of event data as outputted by
#' \code{vast-tools combine}. e.g. each row is an event containing exon
#' metadata, PSI and quality scores values. If \code{expr = TRUE}, then each row
#' is a gene containing two columns of metadata, one column of cRPKM per sample,
#' and, optionally, a second column per sample with read counts.
#' @param expr Set to \code{TRUE} if formatting a cRPKM table. Otherwise, \code{FALSE}.
#' @param counts Set to \code{TRUE} if the cRPKM table has read counts. Otherwise,
#' \code{FALSE}.
#' @param trim_colnames String that must be searched for and trimmed at the end
#' of every sample column in x. Useful to trim the "-cRPKM" suffix from expression
#' tables.
#' @return Data frame with PSI and quality scores. Rownames set as a
#' concatenation of exon metadata delimited by |.
#' @seealso \code{\link{convert_psi}}
#' @export
#' @examples
#' # For example input, see:
#' psi
#' format_table(psi)
#'
#' # For cRPKM
#' format_table(crpkm, expr = TRUE)
format_table <- function(x, expr = FALSE, counts = FALSE, trim_colnames = NULL) {
  if (expr) {
    if (!grepl("^ID$", colnames(x)[1])) {
      stop("Invalid column names. Does your input file contain the correct header?")
    }
    id <- paste(x$NAME,"|",x$ID)
    if (counts) {
      r <- x[,seq(3,ncol(x),by=2)]
    } else {
      r <- x[,3:ncol(x)]
    }
    rownames(r) <- id

    if (!is.null(trim_colnames)){
      str_to_trim <- paste(trim_colnames,"$",sep="")
      colnames(r) <- sapply(colnames(r),function(x) gsub(str_to_trim,"",x))
    }
    return(r)
  }

  if (!grepl("^GENE", colnames(x)[1])) {
    stop("Invalid column names. Does your input file contain the correct header?")
  }
  # Format table to keep only PSIs and convert exon metadata as rownames
  id <- paste(x$COMPLEX, x$GENE, x$COORD, x$LENGTH, sep="|")

  # Extract PSIs
  r <- convert_psi(x[,7:ncol(x)])
  rownames(r) <- id
  return(r)
}
