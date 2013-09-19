# PSI Plotter

source("preprocess_sample_colors.R")

#### Arguments ####
# - Input file
# - Tissue group or Species

print_help <- function() {
  updated <- "2013-Sept-19"
  text <- "**** PSI Plotter ****\nUpdated: 2013-Sept-19\n\nUsage: R CMD BATCH PSI_Plotter.R PSI_Input.tab Tissue_Groups.txt

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

file <- args[1]

# Figure out what tissue data to use 
tissueFile <- args[2]
# if (grepl("mouse|mmu", file, ignore.case=T)) {
#   tissueFile <- "test_data/Tissues.Mmu.txt"
# } else if (grepl("human|hsa", file, ignore.case=T)) {
#   tissueFile <- "../../HONG/Tissues.Hsa.txt"
# } else {
#   stop("Invalid tissue group option - please double check")
# }

if (!file.exists(file))
  stop("Input PSI file doesn't exist!")
if (!file.exists(tissueFile))
  stop("Tissue Group file doesn't exist!")

cat("\n// Input file:", file, "\n")
cat("// Tissue Group file:", tissueFile, "\n")
#########################################################################

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

# assign list of colors
supercolors <- reordered.PSI$col

# Set output file
cat("// Plotting...\n")
outfile <- sub("\\.[^.]*(\\.gz)?$", ".pdf", file)

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

  for (j in 1:ncol(PSIs)){
    abline(v=j,col="grey",lwd=0.3,lty=2)
  }
}
dev.off()

cat("// Done!\n")
cat("// Plots are saved in:", outfile, "\n")
####
