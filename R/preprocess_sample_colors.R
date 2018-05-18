#' Preprocess PSI or cRPKM data frame using configuration file
#'
#' \code{preprocess_sample_colors} re-orders PSI or cRPKM sample columns to a
#' specified order, and it defines sample pools and colors. Order and colors are
#' taken from a pre-defined tab-delimited file (see \code{details}).
#'
#' @details
#' \code{preprocess_sample_colors} depends on a pre-defined "sample inventory"
#' database file in tab-delimited format. This file is species-specific and
#' consists of five columns: \emph{Order}, \emph{SampleName}, \emph{SubroupName}
#' (optional), \emph{GroupName}, \emph{RColorCode}. The header is required in
#' the file. Order of the columns is flexible.
#'
#' For example:
#'
#' \preformatted{Order    SampleName    SubgroupName    GroupName    RColorCode
#' 1        Oocyte_a      Oocyte          EarlyDev     firebrick4
#' 2        Oocyte_b      Oocyte          EarlyDev     firebrick4
#' 3        Embr_4C_a     Embr_4C         EarlyDev     firebrick4
#' 4        Embr_4C_b     Embr_4C         EarlyDev     firebrick4
#' 5        ESC_CGR8      ESC             ESC          firebrick
#' etc..}
#'
#' where:
#' \itemize{
#'  \item{Order: The ordering of the samples from left to right.}
#'  \item{SampleName: Name of the sample. MUST match sample name in input table.}
#'  \item{SubgroupName: Use to define sample pools that will be plotted in the
#'  same data point (see \code{plot_event}, \code{plot_expr} and
#'  \code{plot_multievent}).}
#'  \item{GroupName: Use for plotting the average PSI of samples belonging
#' to the same group. Averages will be calculated from the individual samples,
#' not from the subgroups (to avoid overrepresentation of subgroups with fewer
#' samples).}
#'  \item {RColorCode: Color name as specified by \code{\link{colors}} or hex
#'  color code (\code{#RRGGBB}).}
#' }
#'
#' The \emph{SampleName} must match the column names in \code{data}. It is possible
#' for \code{config} to contain more samples than the \code{data}. In this case,
#' the extra samples will be ignored. It is also possible that
#' \code{config} contains only a subset of the samples in \code{data}. In this
#' case, only the samples specified in the \code{config} will be plotted and everything
#' else is ignored.
#'
#' If a \emph{SampleName} or \emph{SubgroupName} is matched to multiple groups,
#' only the first match will be used. Similarly, if a \emph{GroupName} is matched
#' to multiple \emph{RColorCodes}, the first one will be applied to all elements
#' in the group.
#'
#' To use the \emph{SubgroupNames} in \code{config}, \code{subg} must be set to
#' \code{TRUE} AND \code{config} must contain a \emph{SubgroupName} column. If any
#' of these conditions is not met, one subgroup will be created for each sample,
#' preserving their name and order, and overriding any subgroups in \code{config}.
#'
#' The colors in \code{config} can be overridden by specifiying \code{col}. This
#' was mainly added to support the \code{col} option provided by
#' \code{\link{plot_event}} -- particularly when \code{config} is not provided.
#'
#' This function is also used for formatting cRPKM input data by setting
#' \code{expr = TRUE}.
#'
#' @param data A \emph{n} x \emph{2*m+1} data frame of PSI and quality score values
#' where \emph{n} is the number of AS events and \emph{m} is the number of samples.
#' If \code{expr=TRUE}, a \emph{n} x \emph{m+1} data frame of cRPKM. In both cases,
#' the first column corresponds to the row metadata. Metadata column values must
#' be unique, duplicated values will be discarded.
#' @param config Filename of the configuration file for \code{data}. Also
#' accepts \emph{m} x \emph{4} data frame of the configuration file, or an \emph{m}
#' x \emph{5} data frame if the SubgroupName column is included.
#' @param subg Set to \code{TRUE} to define a subgroup structure using the
#' \code{SubgroupName} column in \code{config}. If \code{FALSE}, or if the file
#' does not contain this column, samples will not be subgrouped (a separate
#' subgroup will be defined for each sample, preserving the sample names). If
#' \code{FALSE} and the file contains a \code{SubgroupName} column, that column
#' will be ignored.
#' @param expr Set to \code{TRUE} if formatting a cRPKM table. Otherwise, \code{FALSE}.
#' @param col Vector of colors with length matching the number of samples. If
#' specified, this will override the color settings specified in \code{config}.
#' @param multi_col Vector of colors with length matching the number of rows in
#' \code{data}. If specified, this can be used to define the color corresponding
#' to each event in \code{plot_multievent()}
#' @return
#' A list containing:
#' \describe{
#'  \item{data}{data frame of PSI/cRPKM values with columns re-ordered}
#'  \item{qual}{data frame of quality scores with columns re-ordered. \code{NULL}
#'  if \code{expr = TRUE}}
#'  \item{sample_order}{data frame with the order corresponding to each sample
#'  name}
#'  \item{subgroup}{data frame with the subgroup corresponding to each sample. If
#'  \code{subg=FALSE}, or \emph{SubgroupName} is not present in \code{config},
#'  a subgroup is made for each sample, preserving sample names.}
#'  \item{subgroup_order}{data frame with the order corresponding to each subgroup.}
#'  \item{group}{data frame with the group corresponding to each subgroup}
#'  \item{group_order}{data frame with the order and color corresponding to
#'  each group.}
#'  \item{multi_col}{if \code{multi_col} was specified, a data frame with the
#'  color coresponding to each event/gene ID. Else, \code{NULL}.}
#'  \item{config}{if a \code{config} was supplied, config data frame summarising
#'  the order-sample-subgroup-group-color relationships described in
#'  \code{sample_order}, \code{subgroup}, \code{subgroup_order}, \code{group}
#'  and \code{group_order}, after correcting for ambiguous relationships. If
#'  \code{col} was supplied, colours in config are overridden with \code{col}. If
#'  no \code{config} was supplied, one will be composed with default parameters.}
#'  \item{original_config}{data frame with the \code{config} supplied to the function}
#'  }
#' @seealso \code{\link{plot_event}}, \code{\link{plot_expr}},
#' \code{\link{plot_multievent}}
#' @export
#' @import dplyr
#' @import readr
#' @import magrittr
#' @examples
#' reorderedpsi <- preprocess_sample_colors(psi, config = config)
#'
#' reorderedcrpkm <- preprocess_sample_colors(crpkm, config = config, expr = TRUE)
preprocess_sample_colors <- function(data,
                                     config,
                                     subg = TRUE,
                                     expr = FALSE,
                                     col = NULL,
                                     multi_col = NULL) {

  R <- list()
  N <- ifelse(expr, ncol(data) - 1, (ncol(data) - 1) / 2)

  if (!is.null(col) && length(col) != N) {
    stop("The length of col does not match the number of samples")
  }

  if(!is.null(multi_col) && length(multi_col) != nrow(data)){
    stop("The length of multi_col does not match the number of events")
  }

  #Filter out non-unique events
  if(length(unique(data$ID)) < length(data$ID)){
    warning("More than one row with the same ID detected. Plotting only one row per ID")
    multi_col <- multi_col[!duplicated(data$ID)]
    data %<>% filter(!duplicated(data$ID))
  }

  #Defining event colours
  #If in multi_event mode and no event colours provided, default ggplot2 will be used
  if(!is.null(multi_col)){
    multi_col <- tibble(ID=data$ID,
                        EventRColorCode=multi_col)
  } else{
    multi_col <- NULL
  }

  #Try to use config from argument
  if(!is.null(config)){

    if (is.character(config)) {
      config <- read_tsv(config)
    }

    original_config <- config

    # check input file [changed to allow flexible config column order, and subgroups]
    if (!all(c("Order", "SampleName", "GroupName", "RColorCode") %in% colnames(config))) {
      stop("Incorrect formatting of headers in config")
    }

    #To use the subgroups in the config, the column must be there and subg==T
    #Otherwise, subgroups equal the individual samples.

    subg <- all("SubgroupName" %in% colnames(config),subg==TRUE)

  } else { #Make our own config and set subg to F

    if (is.null(col)) { #If no colours, use default
      col <- rep("black", N)
    }
    if (expr) {
      data.new <- data[,-1]
      qual.new <- NULL
    } else {
      data.new <- data[, seq(2, ncol(data), 2)]
      qual.new <- data[, seq(3, ncol(data), 2)]
      colnames(qual.new) <- colnames(data.new)
    }

    config <- data.frame("Order"=seq(1,N),
                         "SampleName"=names(data.new),
                         "SubgroupName"=names(data.new),
                         "GroupName" = names(data.new),
                         "RColorCode"=col,
                         stringsAsFactors = F)

    subg <- FALSE
    original_config <- NULL

  }

  #Override colors in config if col argument is there

  if (!is.null(col)) {
    config$RColorCode <- col
  }

  # keep only tissue groups that are present in input data
  # (to take into account samples that might have been excluded)

  config <- config[config$SampleName %in% colnames(data)[-1],]

  if (nrow(config) == 0) {
    stop("No matching samples found in config. Are you using the correct config?")
  }

  # Re-order the PSI table
  config <- config %>%
    arrange(Order) %>%
    mutate(Order=seq(1,nrow(config)))

  new.column.idx <- sapply(config$SampleName,
                           function(x) which(colnames(data) == x))

  data.new <- data[,new.column.idx]

  #If subgroups are not enabled or they are not present in config,
  #then make subgroups equal to samples.

  if(subg==FALSE){
    config <- config %>%
      mutate(SubgroupName=SampleName)
  }

  sample_order <- config %>%
    dplyr::select(Order,SampleName) %>%
    dplyr::rename(SampleOrder="Order") %>%
    as.data.frame()

  subgroup_common <- config %>%
    dplyr::select(Order,SampleName,SubgroupName) %>%
    group_by(SubgroupName) %>%
    dplyr::arrange(Order, .by_group=T)

  subgroup <- subgroup_common %>%
    ungroup() %>%
    dplyr::select(SampleName,SubgroupName) %>%
    as.data.frame()

  subgroup_order <- subgroup_common %>%
    dplyr::select(Order,SubgroupName) %>%
    dplyr::summarise(SubgroupOrder=min(Order)) %>%
    arrange(SubgroupOrder) %>%
    mutate(SubgroupOrder=seq(1:nrow(.))) %>%
    dplyr::select(SubgroupOrder,SubgroupName) %>%
    as.data.frame()

  group_common <- config %>%
    dplyr::select(Order,SubgroupName,GroupName,RColorCode) %>%
    group_by(SubgroupName) %>%
    dplyr::arrange(Order,.by_group=T)

  group <- group_common %>%
    dplyr::select(Order,SubgroupName,GroupName) %>%
    dplyr::summarise(GroupName=first(GroupName)) %>%
    as.data.frame()

  group_order <- group_common %>%
    ungroup() %>%
    group_by(GroupName) %>%
    dplyr::select(Order,GroupName,RColorCode) %>%
    dplyr::summarise(GroupOrder=min(Order),
              RColorCode=first(RColorCode)) %>%
    dplyr::select(GroupOrder,GroupName,RColorCode) %>%
    dplyr::arrange(GroupOrder) %>%
    mutate(GroupOrder=seq(1:nrow(.))) %>%
    as.data.frame()

  config <- sample_order %>%
    left_join(subgroup,by="SampleName") %>%
    left_join(subgroup_order,by="SubgroupName") %>%
    left_join(group,by="SubgroupName") %>%
    left_join(group_order,by="GroupName") %>%
    dplyr::select(Order=SampleOrder,
                 SampleName,
                 SubgroupName,
                 GroupName,
                 RColorCode) %>%
    as.data.frame()

  if (expr) {
    qual.new <- NULL
  } else {
    qual.new <- data[,new.column.idx + 1]
    names(qual.new) <- names(data.new)
  }

  R <- list(data=data.new,
            qual=qual.new,
            sample_order=sample_order,
            subgroup=subgroup,
            subgroup_order=subgroup_order,
            group=group,
            group_order=group_order,
            multi_col=multi_col,
            config=config,
            original_config=original_config)

  return(R)
}
