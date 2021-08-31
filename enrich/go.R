# Title     : GO enrichment for DEGs
# Created by: stl23
# Created on: 2021/8/6
## code from https://blog.csdn.net/sinat_30623997/article/details/79250940?utm_source=blogxgwz1

library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_path <- args[2]
prefix <- args[3]
db <-  args[4]
out_file <- paste0(out_path,'/',prefix,"_GO_out.tsv")
data <- read.table(file=input_file,header=T,sep = "\t")
target_gene_id <- unique(data$GeneID)
if (grepl('^ENSG',target_gene_id[1])){
    target_gene_id_changed <- bitr(target_gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb=db)
    target_gene_id_changed <- target_gene_id_changed$ENTREZID
    }else if (grepl('^NG_',target_gene_id[1])| grepl('^NM_',target_gene_id[1]) |grepl('^NP_',target_gene_id[1]) ){
    target_gene_id_changed <- bitr(target_gene_id, fromType="REFSEQ", toType="ENTREZID", OrgDb=db)
    target_gene_id_changed <- target_gene_id_changed$ENTREZID
    }else if (grepl('^[0-9]+',target_gene_id[1])){
    target_gene_id_changed <- target_gene_id
    }else {
    target_gene_id_changed <- bitr(target_gene_id, fromType="SYMBOL", toType="ENTREZID", OrgDb=db)
    target_gene_id_changed <- target_gene_id_changed$ENTREZID
    }
## out file
ego_ALL <- enrichGO(OrgDb=db,
                   gene = target_gene_id_changed,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   ont = "ALL",
                   readable=TRUE)
write.csv(data.frame(ego_ALL@result),file=out_file,quot=FALSE)
## draw
display_num <- c(15,10,15)
ego_MF <- enrichGO(OrgDb=db,
             gene = target_gene_id_changed,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             qvalueCutoff = 0.05,
             ont = "MF",
             readable=TRUE)
MF_table <- as.data.frame(ego_MF)
if (nrow(MF_table) < display_num[1]){
    ego_result_MF <- MF_table
    }else{
    ego_result_MF <- MF_table[1:display_num[1],]
}
ego_CC <- enrichGO(OrgDb=db,
                   gene = target_gene_id_changed,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
CC_table <- as.data.frame(ego_CC)
if (nrow(CC_table) < display_num[2]){
    ego_result_CC <- CC_table
    }else{
    ego_result_CC <- CC_table[1:display_num[2],]
}
ego_BP <- enrichGO(OrgDb=db,
                   gene = target_gene_id_changed,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
BP_table <- as.data.frame(ego_BP)
if (nrow(BP_table) < display_num[3]){
    ego_result_BP <- BP_table
    }else{
    ego_result_BP <- BP_table[1:display_num[3],]
}
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
type=factor(c(rep("biological process", nrow(ego_result_BP)), rep("cellular component", nrow(ego_result_CC)),
rep("molecular function", nrow(ego_result_MF))), levels=c("biological process", "cellular component", "molecular function"
)))
## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],collapse=" "), "...", sep="")
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

p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")

pdf(file=paste0(out_path,'/',prefix,'_go_enrichment.pdf'))
p
dev.off()

svg(file=paste0(out_path,'/',prefix,'_go_enrichment.svg'))
p
dev.off()

pdf(file=paste0(out_path,'/',prefix,'_go_MF_DAG.pdf'),width = 10,height = 15)
plotGOgraph(ego_MF)
dev.off()
pdf(file=paste0(out_path,'/',prefix,'_go_CC_DAG.pdf'),width = 10,height = 15)
plotGOgraph(ego_CC)
dev.off()
pdf(file=paste0(out_path,'/',prefix,'_go_BP_DAG.pdf'),width = 10,height = 15)
plotGOgraph(ego_BP)
dev.off()
