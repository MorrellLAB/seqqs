# seqqs produces three output files
# 1) nucleotide table, 2) read length by position, 3) a Phred quality matrix by nucleotide position

# nucleotide table plot
# ignoring 'N' and other ambiguities as they are few or absent
pdf(file='~/Desktop/nucleotide_frequency.pdf', 6, 6)
nucl <- read.table('~/Desktop/nucl.txt', header=T)
#nucl <- as.matrix(nucl)
total_bases <- apply(nucl, 1, sum)
A <- nucl$A / total_bases
C <- nucl$C / total_bases
G <- nucl$G / total_bases
T <- nucl$T / total_bases
# legends and figure limits need to be positioned dynamically
plot(A, col='green', ylim=c(0.225, 0.275), cex=0.6, xlab='Nucleotide Position in bp', ylab='Percentage of Nucleotides')
legend(90, 0.275, c('A', 'C', 'G', 'T'), pch=1, col=c('green', 'blue', 'black', 'red'), cex=0.6)
points(C, col='blue', cex=0.6)
points(G, col='black', cex=0.6)
points(T, col='red', cex=0.6)
dev.off()

# need to minimize difference in count because most values are in the 100 bp class
# using log of count
pdf(file='~/Desktop/length.pdf',5,5)
plot(len$pos,log(len$count), xlab='Read Length', ylab='Log of Read Count')
dev.off()

# plot of matrix of quality scores by position
# this needs
pdf(file='~/Desktop/quality.pdf', 8.5, 11)
qual <- read.table('~/Desktop/qual.txt',header=T)
heatmap(as.matrix(qual[apply(qual,2,sum)> 0]),Colv=NA, Rowv=NA, col=gray(9:1/9))
dev.off()