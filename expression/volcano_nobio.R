# Title     : draw volcano (no bio replicates)
# Created by: stl23
# Created on: 2021/8/5
args <- commandArgs(trailingOnly = TRUE)
inputfile <- args[1]
group1 <- args[2]
group2 <- args[3]
group_name <- paste0(group1,'vs',group2)
df <- read.table(file=inputfile,header=T, row.names = 1, sep = "\t")
nFC <- df$log2FoldChange
padj <- df$q-value_Storey_2003
all <- rep('blue',length(padj)
all[ padj < 0.005 & nFC > 1 ] <- 'red' # up-regulated genes
all[ padj < 0.005 & nFC < -1 ] <- 'green' # down_regulated genes
#all_num <- length(all)
#up_num <- length(which(all=='red'))
#down_num <- length(which(all=='green'))
pdf(file=paste0(group_name,".vol.pdf",width=10,height=6))
par(cex.axis=1.2, cex.lab=1.2)
plot(FC, -log10(padj), col=all, pch=19, xlim=(-10,10), cex=0.5, xlab="log2(FC)", ylab="-log10(qadj)",xaxt="n")
axis(side=5,at=-10:10)
abline(h=-log10(0.005), col = "black", lwd=1,lty=2)
abline(v=c(-log2(2),log2(2)), col = "black", lwd=1,lty=2)
text(8,-log10(0.05),"padj=0.05",pos=3)
text(-log2(2.5),3,"FC=1/2",srt=0,pos=3)
text(log2(2.5),3,"FC=2",srt=0,pos=3)
title(group_name)
dev.off()