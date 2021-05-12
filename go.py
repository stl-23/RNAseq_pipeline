import os
import sys

def makeGO(file,group1,group2,out_dir,prefix,OrgDb='org.Hs.eg.db',):
    ## code from https://blog.csdn.net/sinat_30623997/article/details/79250940?utm_source=blogxgwz1
    out = os.path.join(out_dir, prefix + "_GO_out.tsv")
    cmd = """library(clusterProfiler)
data <- read.table(file="%(file)s",header=T,sep = "\t")
target_gene_id <- unique(data$GeneID)
## out file
ego_ALL <- enrichGO(OrgDb="%(OrgDb)s",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "ALL",
                   readable=TRUE)
write.tsv(data.frame(ego_ALL@result),'%(out)s',quote=FALSE)
## draw
display_num <- c(15,10,15)
ego_MF <- enrichGO(OrgDb="%(OrgDb)s",
             gene = target_gene_id,
             qvalueCutoff = 0.05,
             ont = "MF",
             readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_num[1],]
ego_CC <- enrichGO(OrgDb="%(OrgDb)s",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_num[2],]
ego_BP <- enrichGO(OrgDb="%(OrgDb)s",
                   gene = target_gene_id,
                   qvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- as.data.frame(ego_BP)[1:display_num[3],]

go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process"
)))
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}

labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")

p

pdf("%(out_dir)s/%(prefix)s_go_enrichment.pdf")
p
dev.off()

svg("%(out_dir)s/%(prefix)s_go_enrichment.svg")
p
dev.off()
""" %{'file':file,'OrgDb':OrgDb,'out':out,'out_dir':out_dir,'prefix':prefix}

    return cmd