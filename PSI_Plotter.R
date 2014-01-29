# PSI Plotter (aka "Nuno" plots)

source("preprocess_sample_colors.R")

#### Arguments ####
# - Input file
# - Tissue group or Species

print_help <- function() {
  text <- "**** PSI Plotter ****\nUpdated: 2013-Oct-3\n
Usage: Rscript PSI_Plotter.R --args PSI_Input.tab[.gz] Tissue_Groups.txt

Arguments:
    1) Input PSI data - one AS event per row - using the standard PSI format
      e.g. GENE  EVENT  COORD  LENGTH FullCO  COMPLEX  Tissue1_PSI Tissue1_Q ... 
    2) Tissue group file or species (currently supports Hsa or Mmu)
"
  writeLines(text)
}

args <- commandArgs(TRUE)

if (length(args) < 2) {
  print_help()
  stop("Missing arguments")
}

file <- args[2]
tissueFile <- args[3]

if (!file.exists(file))
  stop("Input PSI file doesn't exist!")
if (!file.exists(tissueFile))
  stop("Tissue Group file doesn't exist!")

cat("\n// Input file:", file, "\n")
cat("// Tissue Group file:", tissueFile, "\n")

#### Format input data ####

all_events <- read.csv(file, sep="\t")

format_table <- function(m) {
  id <- paste(m$COMPLEX, m$GENE, m$COORD, m$LENGTH, sep="=")
  psi <- m[,seq(7, ncol(m), 2)]
  rownames(psi) <- id
  return(psi)
}

if (!grepl("^GENE", colnames(all_events)[1])) {
  stop("Invalid column names. Does your input file contain the correct header?")
}
cat("// Brewing some coffee...\n")

cat("// Formatting input data for plotting...\n")
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

cat("//", ncol(PSIs), "samples detected\n")

#### Prepare plotting ####
cat("// Plotting...\n")

# assign list of colors
supercolors <- reordered.PSI$col

# Set output file
outfile <- sub("\\.[^.]*(\\.gz)?$", ".PSI_plots.pdf", file)

par(mfrow=c(1,1),las=2) #3 graphs per row; 2=label always perpendicular to the axis
pdf(outfile,width=8.5,height=5.5)
for(i in 1:nrow(PSIs)){
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

cat("// Done!\n")
cat("//", nrow(PSIs), "plots are saved in:", outfile, "\n")
####
