import os
import sys

def makeKEGG(file,out_dir,prefix,organism='hsa'):
    out = os.path.join(out_dir, prefix + "_KEGG_out.tsv")
    cmd = """library(clusterProfiler)
library(KEGG.db) ## kegg.db was built according to https://github.com/YuLab-SMU/createKEGGdb
data <- read.table(file="%(file)s",header=T,sep = "\t")
target_gene_id <- unique(data$GeneID)
kegg_df <- enrichKEGG(target_gene_id,
                      organism=organism,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      )
write.tsv(data.frame(kegg_df@result),'%(out)s',quot=FALSE)
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
"""%{'file':file,'out':out,}