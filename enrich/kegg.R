# Title     : KEGG enrichment for DEGs
# Created by: stl23
# Created on: 2021/8/6

library(clusterProfiler)
library(KEGG.db) ## kegg.db was built according to https://github.com/YuLab-SMU/createKEGGdb
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_path <- args[2]
prefix <- args[3]
Org <-  args[4]
outfile <- paste0(out_path,'/',prefix,'_KEGG_out.tsv')
data <- read.table(file=input_file,header=T,sep = "\t")
target_gene_id <- unique(data$GeneID)
kegg_df <- enrichKEGG(target_gene_id,
                      organism=Org,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      )
write.tsv(data.frame(kegg_df@result),file=,quot=FALSE)
## draw
display_num <- 20
num <- nrow(as.data.frame(kegg_df))
if (num < display_num){
    barplot(kegg_df, showCategory = num)
    dotplot(kegg_df, showCategory = num)
    dev.off()
    }else{
    barplot(kegg_df, showCategory = display_num)
    dotplot(kegg_df, showCategory = display_num)
    dev.off()
}
