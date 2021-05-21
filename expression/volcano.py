
def volcano(file,group1,group2,biorepeat=True):
    ## bio replicates: padj < 0.05 && log2FC >0 (up regulated); padj < 0.05 && log2FC < 0 (down regulated)
    ## no bio replicates: padj < 0.005 && log2FC >1 (up regulated); padj < 0.005 && log2FC < -1 (down regulated)
    # code from https://www.yunbios.net/R-language-volcano-mapping.html
    if biorepeat:
        cmd = f"""## draw volcano (bio replicates)
df <- read.table(file="{file}",header=T, row.names = 1, sep = "\t")
#p = ggplot(df, aes(log2FoldChange, -1*log10(padj))) 
#p+ geom_point()
nFC <- df$log2FoldChange
padj <- df$padj
all <- rep('blue',length(padj) 
all[ padj < 0.05 & nFC > 0 ] <- 'red' # up-regulated genes
all[ padj < 0.05 & nFC < 0 ] <- 'green' # down_regulated genes
#all_num <- length(all)
#up_num <- length(which(all=='red'))
#down_num <- length(which(all=='green'))
par(cex.axis=1.2, cex.lab=1.2)
plot(nFC, -log10(padj), col=all, pch=19, xlim=(-10,10), cex=0.5, xlab="log2(FC)", ylab="-log10(qadj)",xaxt="n")
axis(side=5,at=-10:10)
abline(h=-log10(0.05), col = "black", lwd=2,lty=2)
text(8,-log10(0.05),"padj=0.05",pos=3)
title("{group1}vs{group2}")
""".format(**locals())
    else:
        cmd = f"""## draw volcano (no bio replicates)
df <- read.table(file="{file}",header=T, row.names = 1, sep = "\t")
nFC <- df$log2FoldChange
padj <- df$q-value_Storey_2003
all <- rep('blue',length(padj) 
all[ padj < 0.005 & nFC > 1 ] <- 'red' # up-regulated genes
all[ padj < 0.005 & nFC < -1 ] <- 'green' # down_regulated genes
#all_num <- length(all)
#up_num <- length(which(all=='red'))
#down_num <- length(which(all=='green'))
par(cex.axis=1.2, cex.lab=1.2)
plot(FC, -log10(padj), col=all, pch=19, xlim=(-10,10), cex=0.5, xlab="log2(FC)", ylab="-log10(qadj)",xaxt="n")
axis(side=5,at=-10:10)
abline(h=-log10(0.005), col = "black", lwd=1,lty=2)
abline(v=c(-log2(2),log2(2)), col = "black", lwd=1,lty=2)
text(8,-log10(0.05),"padj=0.05",pos=3)
text(-log2(2.5),3,"FC=1/2",srt=0,pos=3)
text(log2(2.5),3,"FC=2",srt=0,pos=3)
title("{group1}vs{group2}")
""".format(**locals())

    return cmd