#' Preprocess PSI data frame using configuration file
#'
#' \code{preprocess_sample_colors} re-orders PSI sample columns to a specified order
#' and generates a corresponding color code sequence. The sample order and
#' corresponding colors are taken from a pre-defined tab-delimited file
#' (see \code{details}).
#'
#'@usage preprocess_sample_colors(psi, config)
#'
#' @details
#' \code{preprocess_sample_colors} depends on a pre-defined "sample inventory"
#' database file in tab-delimited format. This file is species-specific and
#' consists of four columns: \emph{Order, SampleName, GroupName, RColorCode}. The
#' header is required in the file.
#'
#' For example:
#'
#' \preformatted{Order    SampleName    GroupName    RColorCode
#' 1        Ooctye        EarlyDev     blue
#' 2        Embr_2C       EarlyDev     red
#' etc..}
#'
#' where:
#' \itemize{
#'  \item Order: The ordering of the samples from left to right.
#'  \item SampleName: Name of the sample. MUST match sample name in input table.
#'  \item Group name: Use for plotting the average PSI of samples belonging
#' to the same group (need to use option -u/--group-means)
#'  \item RColorCode: Color name as specified by \code{colors()} or hex color code (\code{#RRGGBB})
#' }
#'
#' The \emph{SampleName} must match the column names in \code{psi}. It is possible
#' for \code{config} to contain more samples than the \code{psi}. In this case,
#' the extra samples will be ignored. It is also possible that
#' \code{config} contains only a subset of the samples in \code{psi}. In this
#' case, only the samples specified in the \code{config} will be plotted and everything
#' else is ignored.
#'
#' @param data a \emph{n} x \emph{2*m} data frame of PSI and quality score values where \emph{n} is the number of AS events
#' and \emph{m} is the number of samples.
#' @param config filename of the configuration file for \code{data}. Also
#' accepts \emph{m'} x \emph{4} data frame of the configuration file, \code{m <=
#' m}
#' @export
#' @return
#' A list containing:
#'  \item{data}{data frame of PSI values with columns re-ordered}
#'  \item{qual}{data frame of quality scores with columns re-ordered}
#'  \item{col}{vector of colors corresponding to the re-ordered columns}
#'  \item{group.index}{list of column indices corresponding to each \emph{GroupName}}
#'  \item{group.col}{vector of colors corresponding to each \emph{GroupName}}
#' @examples
#' reorderedpsi <- preprocess_sample_colors(mm.psi, config = mm.psi.config)
preprocess_sample_colors <- function(data, config) {
   R <- list()

   if (is.null(config)) {
     mycols <- rep("black", ncol(data) / 2)
     R <- list(data=data[, seq(1, ncol(data), 2)],
               qual=data[, seq(2, ncol(data), 2)],
               col=mycols, group.index=NULL, group.col=NULL)
   } else {
     if (is.character(config)) {
       config <- read.table(config, header = T, sep="\t", comment.char="", 
                            stringsAsFactors=FALSE)
     }

     # check input file
     if (!all(colnames(config) == c("Order", "SampleName", "GroupName", "RColorCode"))) {
       stop("Incorrect formatting of headers in config")
     }

     # check if all samples in input data is in the config
     #unk.samples <- colnames(data) %in% config$SampleName
     #if (!all(unk.samples)) {
     #s <- colnames(data)[!unk.samples]
     #stop(paste("The following samples are not in the tissues config:",
     #paste(s, collapse=", ")))
     #}

     # keep only tissue groups that are present in input data
     # (to take into account samples that might have been excluded)
     config <- config[config$SampleName %in% colnames(data),]

     # Re-order the PSI table
     config <- config[order(config$Order),]
     config$Order <- 1:nrow(config)
     new.column.idx <- sapply(config$SampleName,
                              function(x) which(colnames(data) == x))

     data.new <- data[,new.column.idx]

     # Generate a corresponding color code sequence
     mycols <- config$RColorCode
     names(mycols) <- config$SampleName

     # Store indices of columns for each group
     groups <- unique(config$GroupName)
     mygroups <- list()
     mygroupcol <- rep(NA, length(groups))
     for (i in 1:length(groups)) {
       mygroups[[i]] <- which(colnames(data.new) %in%
                                config[config$GroupName == groups[i],"SampleName"])
       mygroupcol[i] <- as.character(config[config$GroupName == groups[i],
                                        "RColorCode"][1])
     }
     names(mygroups) <- groups
     names(mygroupcol) <- groups

     qual.new <- data[,new.column.idx + 1]
     R <- list(data=data.new,
               qual=qual.new,
               col=mycols,
               group.index=mygroups, group.col=mygroupcol)
   }

   return(R)
}
