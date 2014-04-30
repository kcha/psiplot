#!/usr/bin/env Rscript
#
# PSI Plotter script (aka "Nuno" plots)

source("preprocess_sample_colors.R")

version <- function() {
  return("0.3")
}

print_help <- function() {
  text <- "**** PSI Plotter ****
Script for generating PSI plots across samples (aka Nuno plots).

Usage: ./PSI_Plotter.R PSI_Input.tab[.gz] Tissue_Groups.txt

Arguments:
  1) Input PSI data - one AS event per row - using the standard PSI format
    e.g. GENE  EVENT  COORD  LENGTH FullCO  COMPLEX  Tissue1_PSI Tissue1_Q ... 
  2) Tissue group file or species (currently supports Hsa or Mmu)

Options:
  --version     Display the version    
  --help        This help message

Output:
  A PDF file will be created with one PSI plot per page.

Test run:
  ./PSI_Plotter.R test_data/INCLUSION_LEVELS-ALL3m-Mmu89-SELECTED.test.tab \\
    test_data/Tissues.Mmu.txt
"
  writeLines(text, stderr())
  write(paste("Version:", version()), stderr())
  write("Updated: 2014-04-09", stderr())
}

#### Arguments #################################################################
# - Input file
# - Tissue group or Species

args <- commandArgs(TRUE)

if (args[1] == "--help") {
  print_help()
  stop()
}
if (args[1] == "--version") {
  write(paste("Version:", version()), stderr())
  stop()
}
if (length(args) < 2) {
  print_help()
  stop("Missing arguments")
}

file <- args[length(args)-1]
tissueFile <- args[length(args)]

if (!file.exists(file))
  stop(paste("Input PSI file", file, "doesn't exist!"))
if (!file.exists(tissueFile))
  stop(paste("Tissue Group file", tissueFile, "doesn't exist!"))

write(paste("\n// Input file:", file), stderr())
write(paste("// Tissue Group file:", tissueFile), stderr())

#### Format input data #########################################################

all_events <- read.csv(file, sep="\t")

convert_psi <- function(t) {
  # Helper function to filter and return PSI values
  # PSIs are converted to NA if first coverage code is 'N'
  # e.g. PSI=100, Coverage=N,N,N,OK,S ---> PSI=NA
  #
  # Input: original PSI plus quality scores WITHOUT the first 7 columns

  stopifnot(ncol(t) %% 2 == 0)
  psi <- t
  
  for (i in seq(1, ncol(psi), 2)) {
    cov <- strsplit(as.character(psi[,i+1]), split = ",")
    cov <- sapply(cov, "[", 1)
    
    na <- which(cov == "N")
    if (length(na) > 0) 
      psi[na, i] <- NA
  }
  return(psi[, seq(1, ncol(psi), 2)])
}

format_table <- function(m) {
  # Format table to keep only PSIs and convert exon metadata as rownames
  id <- paste(m$COMPLEX, m$GENE, m$COORD, m$LENGTH, sep="=")
  
  # Extract PSIs
  psi <- convert_psi(m[,7:ncol(m)])
  rownames(psi) <- id
  return(psi)
}

if (!grepl("^GENE", colnames(all_events)[1])) {
  stop("Invalid column names. Does your input file contain the correct header?")
}
write("// Brewing some coffee...", stderr())

write("// Formatting input data for plotting...", stderr())
PSIs <- format_table(all_events)
# Call function to re-order columns of PSI data
#
# returns a list containing four elements:
#   data        - the PSI data with sample columsn re-ordered
#   col         - vector of colours that will be plotted
#   group.index - list of indices for each sample group (e.g. ESC, Neural, etc.)
#   group.col   - corresponding color for sample group
reordered.PSI <- preprocess_sample_colors(PSIs, tissueFile)
PSIs <- as.matrix(reordered.PSI$data)
ALLev <- row.names(PSIs)
samples <- colnames(PSIs)

write(paste("//", ncol(PSIs), "samples detected"), stderr())

#### Prepare plotting ##########################################################
write("// Plotting...", stderr())

# assign list of colors
supercolors <- reordered.PSI$col

# Set output file
outfile <- sub("\\.[^.]*(\\.gz)?$", ".PSI_plots.pdf", file)

pdf(outfile, width = 8.5, height = 5.5)
par(mfrow = c(1,1), las = 2) #3 graphs per row; 2=label always perpendicular to the axis
for (i in 1:nrow(PSIs)) {
  plot(as.numeric(PSIs[i,]),
       col=supercolors,
       pch=20,
       main=rownames(PSIs)[i],
       ylab="PSI", xlab="", xaxt="n",
       ylim=c(1,100),
       cex=0.8, cex.main=0.9, cex.axis=0.8)
  axis(1, at=seq(1, ncol(PSIs), by=1), labels = FALSE)
  text(seq(1, ncol(PSIs), by=1), 
       par("usr")[3] - 3.5, 
       labels = samples, 
       srt = 45, adj=c(1,1), xpd = TRUE,cex=0.5)
  
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["ESC"]] ], na.rm=TRUE), 
         col=reordered.PSI$group.col["ESC"], lwd=0.5)
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["Neural"]] ], na.rm=TRUE),
         col=reordered.PSI$group.col["Neural"], lwd=0.5)
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["Muscle"]] ], na.rm=TRUE),
         col=reordered.PSI$group.col["Muscle"], lwd=0.5)
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["Tissues"]] ], na.rm=TRUE),
         col=reordered.PSI$group.col["Tissues"], lwd=0.5)

  abline(v=1:ncol(PSIs), col="grey", lwd=0.3, lty=2)
  abline(h=seq(0,100,10), col="grey", lwd=0.3, lty=2)
}
dev.off()

write("// Done!\n", stderr())
write(paste("//", nrow(PSIs), "plots are saved in:", outfile), stderr())
####
