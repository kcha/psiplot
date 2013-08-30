#All events
all_events=read.table("test_data/INCLUSION_LEVELS-ALL3m-Mmu89-PSI_TEST.tab.gz", header=T, sep="\t",row.names=1)


#### SG's filtering of events --- NOT USED HERE IN THIS TEST EXAMPLE ####
#Paramters
#Infile
# infile="events.txt"
#Outfile
outfile="events_out"

#List of query events
# chosen_events <- scan(infile, what="", sep="\n")

#Indeces of events list in "all events table"
# indeces=grep(paste(chosen_events,collapse="|"),as.character(row.names(all_events)))

#Plot chosen events
# PSIs_table=all_events[indeces,]

#########################################################################

source("preprocess_sample_colors.R")

# call function to re-order columns of PSI data
#
# returns a list containing four elements:
#   data        - the PSI data with sample columsn re-ordered
#   col         - vector of colours that will be plotted
#   group.index - list of indices for each sample group (e.g. ESC, Neural, etc.)
#   group.col   - corresponding color for sample group
reordered.PSI <- preprocess_sample_colors(all_events, "test_data/Tissues.Mmu.txt")
PSIs <- as.matrix(reordered.PSI$data)
ALLev <- row.names(PSIs)
samples <- colnames(PSIs)

# assign list of colors
supercolorsH <- reordered.PSI$col

# cl <- colors()
# supercolorsM=c(rep(cl[32],each=5),rep("red",each=6),rep(cl[506],each=2),rep(cl[503],each=2),
#                rep(cl[375],each=6),rep(cl[655],each=2),rep(cl[145],each=5),rep("gold",each=3),rep(cl[623],each=2),
#                rep(cl[635],each=3),rep(cl[639],each=4),"magenta",cl[639],"magenta",rep("blue",each=8),
#                rep("orange",each=2),"magenta",cl[50],"magenta",rep(cl[50],each=4),rep("green",each=4),
#                rep(cl[600],each=2),rep(cl[599],each=2),rep("black",each=15),rep("orange",each=5))
#earlydev,ESC,EpiSC,iPS
#mid_dev,PGC,Spermatog,cyte,tid,zoa,Sertoli,Testis,EmbrLimb,EmbrLiver
#NPCs,Embr_Brain,nSR100,Neural,
#N2A,nsR100-KD, C2C12,Mbnl,MyoblDif,Muscle,Heart
#TS,Placenta,Tissues,CLines

par(mfrow=c(1,1),las=2) #3 graphs per row; 2=label always perpendicular to the axis
pdf(paste(outfile,".pdf",sep=""),width=8.5,height=5.5)
for(i in 1:nrow(PSIs)){
  plot(as.numeric(PSIs[i,]),col=supercolorsH,pch=20,main=rownames(PSIs)[i],ylab="PSI",xlab="",xaxt="n",
     ylim=c(1,100),cex=0.8,cex.main=0.9,cex.axis=0.8)
  axis(1, at=seq(1, ncol(PSIs), by=1), labels = FALSE)
  text(seq(1, ncol(PSIs), by=1), par("usr")[3] - 3.5, labels = samples, srt = 45, adj=c(1,1), xpd = TRUE,cex=0.5)
  
  # *NEW* way of drawing horizontal lines.
  # Overall, it is the same style as before, just a little bit longer 
  # and complex.
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["ESC"]] ], na.rm=TRUE), 
         col=reordered.PSI$group.col["ESC"], lwd=0.5)
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["Neural"]] ], na.rm=TRUE),
         col=reordered.PSI$group.col["Neural"], lwd=0.5)
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["Muscle"]] ], na.rm=TRUE),
         col=reordered.PSI$group.col["Muscle"], lwd=0.5)
  abline(h=mean(PSIs[i, reordered.PSI$group.index[["Tissues"]] ], na.rm=TRUE),
         col=reordered.PSI$group.col["Tissues"], lwd=0.5)
  
  # The old way (for reference)
#   abline(h=mean(PSIs[i,6:11],na.rm=TRUE),col="red",lwd=0.5)
#   abline(h=mean(PSIs[i,44:51],na.rm=TRUE),col="blue",lwd=0.5)  
#   abline(h=mean(PSIs[i,61:64],na.rm=TRUE),col="green",lwd=0.5)
#   abline(h=mean(PSIs[i,67:88],na.rm=TRUE),col="black",lwd=0.5)
  for (j in 1:ncol(PSIs)){
    abline(v=j,col="grey",lwd=0.3,lty=2)
  }
}
dev.off()

####