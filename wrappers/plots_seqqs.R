# seqqs produces three output files
# 1) nucleotide table, 2) read length by position, 3) a Phred quality matrix by nucleotide position

#   To take command line arguments
args <- commandArgs(TRUE)
#   This creates a vector of character strings for arguments
#   we will just take two arguments here, the stats directory
#   and the sample name
statsdir <- args[1]
samplename <- args[2]

# nucleotide table plot
#   print out some messages
print(paste("Creating raw nucleotide composition plot for ", samplename, " Forward Reads..."), sep="")
#   Build the name of the file that will be read here
raw.forward.nucl.name <- paste(statsdir, "/raw_", samplename, "_R1_nucl.txt", sep="")
#   and build the output name
raw.forward.nucl.output <- paste(statsdir, "/plots/Raw_", samplename, "_R1_NuclPlot.pdf", sep="")
# ignoring 'N' and other ambiguities as they are few or absent
pdf(file=raw.forward.nucl.output, 6, 6)
nucl <- read.table(raw.forward.nucl.name, header=TRUE)
#   Sum across rows to get the total number of calls at each position
total_bases <- apply(nucl, 1, sum)
#   Calculate the proportions of each nucleotide
#   IDEALLY these are around 25% each
A <- nucl$A / total_bases
C <- nucl$C / total_bases
G <- nucl$G / total_bases
T <- nucl$T / total_bases
#   What is the maximum proportion? Add 0.02 to it for the upper bound of the plot
upper_bound <- max(A, C, G, T) + 0.02
#   Same with the lower bound, but subtract 0.02
lower_bound <- min(A, C, G, T) - 0.02
# legends and figure limits need to be positioned dynamically
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.6, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main=paste("Base Composition in Raw ", samplename, " Forward Reads", sep=""))
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=1, col=c('green', 'blue', 'black', 'red'), cex=0.6)
points(C, col='blue', cex=0.6)
points(G, col='black', cex=0.6)
points(T, col='red', cex=0.6)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=2)
dev.off()

# need to minimize difference in count because most values are in the 100 bp class
# using log of count
print(paste("Creating raw read length distribution plot for ", samplename, " Forward Reads..."), sep="")
raw.forward.len.name <- paste(statsdir, "/raw_", samplename, "_R1_len.txt", sep="")
raw.forward.len.output <- paste(statsdir, "/plots/Raw_", samplename, "_R1_LengthDist.pdf", sep="")
pdf(file=raw.forward.len.output,6, 6)
len <- read.table(raw.forward.len.name, header=TRUE)
plot(len$pos,log(len$count), xlab='Read Length', ylab='Log of Read Count', main=paste("Distribution of Read Length in Raw ", samplename, " Forward Reads", sep=""))
dev.off()

# plot of matrix of quality scores by position
# this needs
print(paste("Creating raw base quality heatmap for ", samplename, " Forward Reads..."), sep="")
raw.forward.qual.name <- paste(statsdir, "/raw_", samplename, "_R1_qual.txt", sep="")
raw.forward.qual.output <- paste(statsdir, "/plots/Raw_", samplename, "_R1_QualityHeatMap.pdf", sep="")
pdf(file=raw.forward.qual.output, 8.5, 11)
qual <- read.table(raw.forward.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Raw ", samplename, " Forward Reads", sep=""))
dev.off()


#####
#   The same for reverse
#####
print(paste("Creating raw nucleotide composition plot for ", samplename, " Reverse Reads..."), sep="")
#   Build the name of the file that will be read here
raw.Reverse.nucl.name <- paste(statsdir, "/raw_", samplename, "_R2_nucl.txt", sep="")
#   and build the output name
raw.Reverse.nucl.output <- paste(statsdir, "/plots/Raw_", samplename, "_R2_NuclPlot.pdf", sep="")
# ignoring 'N' and other ambiguities as they are few or absent
pdf(file=raw.Reverse.nucl.output, 6, 6)
nucl <- read.table(raw.Reverse.nucl.name, header=TRUE)
#   Sum across rows to get the total number of calls at each position
total_bases <- apply(nucl, 1, sum)
#   Calculate the proportions of each nucleotide
#   IDEALLY these are around 25% each
A <- nucl$A / total_bases
C <- nucl$C / total_bases
G <- nucl$G / total_bases
T <- nucl$T / total_bases
#   What is the maximum proportion? Add 0.02 to it for the upper bound of the plot
upper_bound <- max(A, C, G, T) + 0.02
#   Same with the lower bound, but subtract 0.02
lower_bound <- min(A, C, G, T) - 0.02
# legends and figure limits need to be positioned dynamically
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.6, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main=paste("Base Composition in Raw ", samplename, " Reverse Reads", sep=""))
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=1, col=c('green', 'blue', 'black', 'red'), cex=0.6)
points(C, col='blue', cex=0.6)
points(G, col='black', cex=0.6)
points(T, col='red', cex=0.6)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=2)
dev.off()

# need to minimize difference in count because most values are in the 100 bp class
# using log of count
print(paste("Creating raw read length distribution plot for ", samplename, " Reverse Reads..."), sep="")
raw.Reverse.len.name <- paste(statsdir, "/raw_", samplename, "_R2_len.txt", sep="")
raw.Reverse.len.output <- paste(statsdir, "/plots/Raw_", samplename, "_R2_LengthDist.pdf", sep="")
pdf(file=raw.Reverse.len.output,6, 6)
len <- read.table(raw.Reverse.len.name, header=TRUE)
plot(len$pos,log(len$count), xlab='Read Length', ylab='Log of Read Count', main=paste("Distribution of Read Length in Raw ", samplename, " Reverse Reads", sep=""))
dev.off()

# plot of matrix of quality scores by position
# this needs
print(paste("Creating raw base quality heatmap for ", samplename, " Reverse Reads..."), sep="")
raw.Reverse.qual.name <- paste(statsdir, "/raw_", samplename, "_R2_qual.txt", sep="")
raw.Reverse.qual.output <- paste(statsdir, "/plots/Raw_", samplename, "_R2_QualityHeatMap.pdf", sep="")
pdf(file=raw.Reverse.qual.output, 8.5, 11)
qual <- read.table(raw.Reverse.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Raw ", samplename, " Reverse Reads", sep=""))
dev.off()


#####
#   And then for the trimmed reads
#####
print(paste("Creating Trimmed nucleotide composition plot for ", samplename, " Forward Reads..."), sep="")
#   Build the name of the file that will be read here
Trimmed.forward.nucl.name <- paste(statsdir, "/Trimmed_", samplename, "_R1_nucl.txt", sep="")
#   and build the output name
Trimmed.forward.nucl.output <- paste(statsdir, "/plots/Trimmed_", samplename, "_R1_NuclPlot.pdf", sep="")
# ignoring 'N' and other ambiguities as they are few or absent
pdf(file=Trimmed.forward.nucl.output, 6, 6)
nucl <- read.table(Trimmed.forward.nucl.name, header=TRUE)
#   Sum across rows to get the total number of calls at each position
total_bases <- apply(nucl, 1, sum)
#   Calculate the proportions of each nucleotide
#   IDEALLY these are around 25% each
A <- nucl$A / total_bases
C <- nucl$C / total_bases
G <- nucl$G / total_bases
T <- nucl$T / total_bases
#   What is the maximum proportion? Add 0.02 to it for the upper bound of the plot
upper_bound <- max(A, C, G, T) + 0.02
#   Same with the lower bound, but subtract 0.02
lower_bound <- min(A, C, G, T) - 0.02
# legends and figure limits need to be positioned dynamically
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.6, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main=paste("Base Composition in Trimmed ", samplename, " Forward Reads", sep=""))
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=1, col=c('green', 'blue', 'black', 'red'), cex=0.6)
points(C, col='blue', cex=0.6)
points(G, col='black', cex=0.6)
points(T, col='red', cex=0.6)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=2)
dev.off()

# need to minimize difference in count because most values are in the 100 bp class
# using log of count
print(paste("Creating Trimmed read length distribution plot for ", samplename, " Forward Reads..."), sep="")
Trimmed.forward.len.name <- paste(statsdir, "/Trimmed_", samplename, "_R1_len.txt", sep="")
Trimmed.forward.len.output <- paste(statsdir, "/plots/Trimmed_", samplename, "_R1_LengthDist.pdf", sep="")
pdf(file=Trimmed.forward.len.output,6, 6)
len <- read.table(Trimmed.forward.len.name, header=TRUE)
plot(len$pos,log(len$count), xlab='Read Length', ylab='Log of Read Count', main=paste("Distribution of Read Length in Trimmed ", samplename, " Forward Reads", sep=""))
dev.off()

# plot of matrix of quality scores by position
# this needs
print(paste("Creating Trimmed base quality heatmap for ", samplename, " Forward Reads..."), sep="")
Trimmed.forward.qual.name <- paste(statsdir, "/Trimmed_", samplename, "_R1_qual.txt", sep="")
Trimmed.forward.qual.output <- paste(statsdir, "/plots/Trimmed_", samplename, "_R1_QualityHeatMap.pdf", sep="")
pdf(file=Trimmed.forward.qual.output, 8.5, 11)
qual <- read.table(Trimmed.forward.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Trimmed ", samplename, " Forward Reads", sep=""))
dev.off()


#####
#   The same for reverse
#####
print(paste("Creating Trimmed nucleotide composition plot for ", samplename, " Reverse Reads..."), sep="")
#   Build the name of the file that will be read here
Trimmed.Reverse.nucl.name <- paste(statsdir, "/Trimmed_", samplename, "_R2_nucl.txt", sep="")
#   and build the output name
Trimmed.Reverse.nucl.output <- paste(statsdir, "/plots/Trimmed_", samplename, "_R2_NuclPlot.pdf", sep="")
# ignoring 'N' and other ambiguities as they are few or absent
pdf(file=Trimmed.Reverse.nucl.output, 6, 6)
nucl <- read.table(Trimmed.Reverse.nucl.name, header=TRUE)
#   Sum across rows to get the total number of calls at each position
total_bases <- apply(nucl, 1, sum)
#   Calculate the proportions of each nucleotide
#   IDEALLY these are around 25% each
A <- nucl$A / total_bases
C <- nucl$C / total_bases
G <- nucl$G / total_bases
T <- nucl$T / total_bases
#   What is the maximum proportion? Add 0.02 to it for the upper bound of the plot
upper_bound <- max(A, C, G, T) + 0.02
#   Same with the lower bound, but subtract 0.02
lower_bound <- min(A, C, G, T) - 0.02
# legends and figure limits need to be positioned dynamically
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.6, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main=paste("Base Composition in Trimmed ", samplename, " Reverse Reads", sep=""))
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=1, col=c('green', 'blue', 'black', 'red'), cex=0.6)
points(C, col='blue', cex=0.6)
points(G, col='black', cex=0.6)
points(T, col='red', cex=0.6)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=2)
dev.off()

# need to minimize difference in count because most values are in the 100 bp class
# using log of count
print(paste("Creating Trimmed read length distribution plot for ", samplename, " Reverse Reads..."), sep="")
Trimmed.Reverse.len.name <- paste(statsdir, "/Trimmed_", samplename, "_R2_len.txt", sep="")
Trimmed.Reverse.len.output <- paste(statsdir, "/plots/Trimmed_", samplename, "_R2_LengthDist.pdf", sep="")
pdf(file=Trimmed.Reverse.len.output,6, 6)
len <- read.table(Trimmed.Reverse.len.name, header=TRUE)
plot(len$pos,log(len$count), xlab='Read Length', ylab='Log of Read Count', main=paste("Distribution of Read Length in Trimmed ", samplename, " Reverse Reads", sep=""))
dev.off()

# plot of matrix of quality scores by position
# this needs
print(paste("Creating Trimmed base quality heatmap for ", samplename, " Reverse Reads..."), sep="")
Trimmed.Reverse.qual.name <- paste(statsdir, "/Trimmed_", samplename, "_R2_qual.txt", sep="")
Trimmed.Reverse.qual.output <- paste(statsdir, "/plots/Trimmed_", samplename, "_R2_QualityHeatMap.pdf", sep="")
pdf(file=Trimmed.Reverse.qual.output, 8.5, 11)
qual <- read.table(Trimmed.Reverse.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Raw ", samplename, " Reverse Reads", sep=""))
dev.off()
