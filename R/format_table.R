#' Format input event data
#'
#' Convert low/bad quality PSI values and convert event metadata as a first ID
#' column. Prepares input event/GE data for plotting. Calls \code{\link{convert_psi}}.
#'
#' @param x A data frame of event data as outputted by \code{vast-tools combine}.
#' e.g. each row is an event containing exon metadata, PSI and quality scores
#' values. If \code{expr = TRUE}, then each row is a gene containing two columns
#' of metadata, one column of cRPKM per sample, and, optionally, a second column
#' per sample with read counts.
#' @param qual String indicating the minimun \emph{vast-tools} quality score
#' for the PSI to be accepted. Defaults to \code{'VLOW'}. See the
#' \href{https://github.com/vastgroup/vast-tools/blob/master/README.md}{vast-tools
#' documentation} for details.
#' @param expr Set to \code{TRUE} if formatting a cRPKM table. Otherwise,
#' \code{FALSE}.
#' @param counts Set to \code{TRUE} if the cRPKM table has read counts. Otherwise,
#' \code{FALSE}.
#' @param trim_colnames String that must be searched for and trimmed at the end
#' of every sample column in x. Useful to trim the "-cRPKM" suffix from expression
#' tables. If no string must be trimmed, leave as \code{FALSE}.
#' @param short_ids Set to \code{TRUE} to make the metadata column shorter in the
#' output, by including only the event ID (for events) or the gene ID (if run
#' with \code{expr==TRUE}).
#' @return A data frame. If run with \code{expr=FALSE} (default), each row is an
#' event, the first column (\code{ID}) contains a concatenation of the event
#' metadata delimited by |, and there are two more columns per sample, with PSI
#' and quality scores values. If run with (\code{expr=TRUE}), each row is a gene,
#' the \code{ID} column contains the gene metadata, and there is one more column
#' per sample with the cRPKM value.
#' @seealso \code{\link{convert_psi}}
#' @export
#' @import dplyr
#' @importFrom stats setNames
#' @examples
#' # For example input, see:
#' psi
#' format_table(psi)
#' # For cRPKM
#' crpkm
#' format_table(crpkm, expr = TRUE)
#'
#' # For cRPKM with read counts and the "-cRPKM" suffix in sample columns:
#' crpkm_counts
#' format_table(crpkm_counts, expr = TRUE, counts = TRUE, trim_colnames = "-cRPKM")
#'
#' # To keep only event IDs/gene IDs as metadata:
#' psi
#' format_table(psi,short_ids = TRUE)
#'
format_table <- function(x,
                         qual = c("VLOW","N","LOW","OK","SOK"),
                         expr = FALSE,
                         counts = FALSE,
                         trim_colnames = NULL,
                         short_ids = FALSE) {

  if (expr) {
    if (!grepl("^ID$", colnames(x)[1])) {
      stop("Invalid column names. Does your input file contain the correct header?")
    }

    #Use only the ID column as gene ID if short_ids==T. Else, paste all metadata
    if(short_ids==FALSE){
      r <- x %>% mutate(ID=paste(NAME,ID,sep="|"))
    } else{
      r <- x
    }

    if (counts) {
      r <- r[,c(1,seq(3,ncol(x),by=2))]
    } else {
      r <- r[,c(1,3:ncol(x))]
    }

    if (!is.null(trim_colnames)){
      str_to_trim <- paste(trim_colnames,"$",sep="") #Regexp for end of line
      r <- setNames(r,
                    c(colnames(r)[1],gsub(str_to_trim,"",colnames(r)[-1])))

    }
    return(r)
  }

  if (!grepl("^GENE", colnames(x)[1])) {
    stop("Invalid column names. Does your input file contain the correct header?")
  }

  # Format table to keep only PSIs an ID column.
  # If short_ids==T, keep only event ID. Else, paste all metadata as ID

  if(short_ids==T){
    id <- x$EVENT
  } else{
    id <- paste(x$COMPLEX, x$GENE, x$EVENT, x$COORD, x$LENGTH, sep="|")
  }

  # Extract PSIs

  r <- convert_psi(x[,7:ncol(x)],qual=qual)

  r <- r %>% mutate(ID=id) %>%
    select(ID,colnames(r))

  return(r)
}
