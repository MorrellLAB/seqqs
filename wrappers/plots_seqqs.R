# seqqs produces three output files
# 1) nucleotide table, 2) read length by position, 3) a Phred quality matrix by nucleotide position

#   To take command line arguments
args <- commandArgs(TRUE)
#   This creates a vector of character strings for arguments
#   we will just take two arguments here, the stats directory
#   and the sample name
statsdir <- args[1]
samplename <- args[2]

#   Start a plot file for forward
outputfile <- paste(statsdir, "/", samplename, "_Forward_SeqqsPlots.pdf", sep="")
pdf(file=outputfile, width=8.8, 11)

#   Set it up so that we have multiple plots on one graphics device
par(mfrow=c(2, 2))
#   They fill top-to-bottom, left-to-right

#####
#   Base composition plots
#####
#       Forward before trimming
#   Build the name of the file that will be read here
raw.forward.nucl.name <- paste(statsdir, "/raw_", samplename, "_R1_nucl.txt", sep="")
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
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.3, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main="Raw\nBase Composition", pch=19)
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=19, col=c('green', 'blue', 'black', 'red'), cex=0.7)
points(C, col='blue', cex=0.3, pch=19)
points(G, col='black', cex=0.3, pch=19)
points(T, col='red', cex=0.3, pch=19)
lines(A, col="green", type="l", lwd=0.5)
lines(C, col='blue', type="l", lwd=0.5)
lines(G, col='black', type="l", lwd=0.5)
lines(T, col='red', type="l", lwd=0.5)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=1)

#       Forward after trimming
Trimmed.forward.nucl.name <- paste(statsdir, "/Trimmed_", samplename, "_R1_nucl.txt", sep="")
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
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.3, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main="Cleaned\nBase Composition", pch=19)
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=19, col=c('green', 'blue', 'black', 'red'), cex=0.7)
points(C, col='blue', cex=0.3, pch=19)
points(G, col='black', cex=0.3, pch=19)
points(T, col='red', cex=0.3, pch=19)
lines(A, col="green", type="l", lwd=0.5)
lines(C, col='blue', type="l", lwd=0.5)
lines(G, col='black', type="l", lwd=0.5)
lines(T, col='red', type="l", lwd=0.5)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=1)

#####
#   Length distribution
#####
#       Before trimming
raw.forward.len.name <- paste(statsdir, "/raw_", samplename, "_R1_len.txt", sep="")
len <- read.table(raw.forward.len.name, header=TRUE)
counts <- log(len$count)
counts[counts == -Inf] <- 0
#   The step size for the bar plot
step <- 5
#len <- len[len$count != 0,]
at <- barplot(counts, xlab='Read Length', ylab='Log of Read Count', main="Length Distribution", col="blue", space=rep(0, length(len)), border="white")
at.labs <- seq_along(at)
at.labs <- c(at.labs[seq(1, length(at.labs), step)], length(at.labs))
at.pos <- c(at[seq(1, length(at), step)], at[length(at)])
axis(1, at=at.pos, labels=as.character(at.labs))

#       After trimming
Trimmed.forward.len.name <- paste(statsdir, "/Trimmed_", samplename, "_R1_len.txt", sep="")
len <- read.table(Trimmed.forward.len.name, header=TRUE)
counts <- log(len$count)
counts[counts == -Inf] <- 0
#   The step size for the bar plot
step <- 5
#len <- len[len$count != 0,]
at <- barplot(counts, xlab='Read Length', ylab='Log of Read Count', main="Length Distribution", col="blue", space=rep(0, length(len)), border="white")
at.labs <- seq_along(at)
at.labs <- c(at.labs[seq(1, length(at.labs), step)], length(at.labs))
at.pos <- c(at[seq(1, length(at), step)], at[length(at)])
axis(1, at=at.pos, labels=as.character(at.labs))

dev.off()


#####
#   Quality heatmaps
#####
#       Before trimming
raw.forward.qual.name <- paste(statsdir, "/raw_", samplename, "_R1_qual.txt", sep="")
raw.forward.output <- paste(statsdir, "/Raw_", samplename, "_ForwardQuality.pdf", sep="")
pdf(file=raw.forward.output, 8, 8)
qual <- read.table(raw.forward.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Raw ", samplename, " Forward Reads", sep=""))
dev.off()

#       After trimming
Trimmed.forward.qual.name <- paste(statsdir, "/Trimmed_", samplename, "_R1_qual.txt", sep="")
Trimmed.forward.output <- paste(statsdir, "/Trimmed_", samplename, "_ForwardQuality.pdf", sep="")
pdf(file=Trimmed.forward.output, 8, 8)
qual <- read.table(Trimmed.forward.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Trimmed ", samplename, " Forward Reads", sep=""))


#####
#   Reverse plots
#####
outputfile <- paste(statsdir, "/", samplename, "_Reverse_SeqqsPlots.pdf", sep="")
pdf(file=outputfile, width=8.8, height=11)
#   Set it up so that we have multiple plots on one graphics device
par(mfrow=c(2, 2))

#####
#   Base composition plots
#####
#   Before trimming
#   Build the name of the file that will be read here
raw.Reverse.nucl.name <- paste(statsdir, "/raw_", samplename, "_R2_nucl.txt", sep="")
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
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.3, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main="Raw\nBase Composition", pch=19)
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=19, col=c('green', 'blue', 'black', 'red'), cex=0.7)
points(C, col='blue', cex=0.3, pch=19)
points(G, col='black', cex=0.3, pch=19)
points(T, col='red', cex=0.3, pch=19)
lines(A, col="green", type="l", lwd=0.5)
lines(C, col='blue', type="l", lwd=0.5)
lines(G, col='black', type="l", lwd=0.5)
lines(T, col='red', type="l", lwd=0.5)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=1)

#   After trimming
Trimmed.Reverse.nucl.name <- paste(statsdir, "/Trimmed_", samplename, "_R2_nucl.txt", sep="")
#   and build the output name
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
plot(A, col='green', ylim=c(lower_bound, upper_bound), cex=0.3, xlab='Position in Read (bp)', ylab='Percentage of Nucleotides', main="Cleaned\nBase Composition", pch=19)
legend(90, upper_bound, c('A', 'C', 'G', 'T'), pch=19, col=c('green', 'blue', 'black', 'red'), cex=0.7)
points(C, col='blue', cex=0.3, pch=19)
points(G, col='black', cex=0.3, pch=19)
points(T, col='red', cex=0.3, pch=19)
lines(A, col="green", type="l", lwd=0.5)
lines(C, col='blue', type="l", lwd=0.5)
lines(G, col='black', type="l", lwd=0.5)
lines(T, col='red', type="l", lwd=0.5)
#   Add a line at 0.25
abline(h=0.25, col="orange", lwd=1)

#####
#   Length distribution
#####
#       Before trimming
raw.Reverse.len.name <- paste(statsdir, "/raw_", samplename, "_R2_len.txt", sep="")
len <- read.table(raw.Reverse.len.name, header=TRUE)
counts <- log(len$count)
counts[counts == -Inf] <- 0
#   The step size for the bar plot
step <- 5
#len <- len[len$count != 0,]
at <- barplot(counts, xlab='Read Length', ylab='Log of Read Count', main="Length Distribution", col="blue", space=rep(0, length(len)), border="white")
at.labs <- seq_along(at)
at.labs <- c(at.labs[seq(1, length(at.labs), step)], length(at.labs))
at.pos <- c(at[seq(1, length(at), step)], at[length(at)])
axis(1, at=at.pos, labels=as.character(at.labs))

#       After trimming
Trimmed.Reverse.len.name <- paste(statsdir, "/Trimmed_", samplename, "_R2_len.txt", sep="")
len <- read.table(Trimmed.Reverse.len.name, header=TRUE)
counts <- log(len$count)
counts[counts == -Inf] <- 0
#   The step size for the bar plot
step <- 5
#len <- len[len$count != 0,]
at <- barplot(counts, xlab='Read Length', ylab='Log of Read Count', main="Length Distribution", col="blue", space=rep(0, length(len)), border="white")
at.labs <- seq_along(at)
at.labs <- c(at.labs[seq(1, length(at.labs), step)], length(at.labs))
at.pos <- c(at[seq(1, length(at), step)], at[length(at)])
axis(1, at=at.pos, labels=as.character(at.labs))
dev.off()

#####
#   Quality heatmaps
#####

#   Before trimming
raw.Reverse.qual.name <- paste(statsdir, "/raw_", samplename, "_R2_qual.txt", sep="")
raw.Reverse.output <- paste(statsdir, "/Raw_", samplename, "_ReverseQuality.pdf", sep="")
pdf(file=raw.Reverse.output, 8, 8)
qual <- read.table(raw.Reverse.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Raw ", samplename, " Reverse Reads", sep=""))
dev.off()

#   After trimming
Trimmed.Reverse.qual.name <- paste(statsdir, "/Trimmed_", samplename, "_R2_qual.txt", sep="")
Trimmed.Reverse.output <- paste(statsdir, "/Trimmed_", samplename, "_ReverseQuality.pdf", sep="")
pdf(file=Trimmed.Reverse.output, 8, 8)
qual <- read.table(Trimmed.Reverse.qual.name,header=TRUE)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9), main=paste("Quality Heatmap in Trimmed ", samplename, " Reverse Reads", sep=""))
dev.off()

