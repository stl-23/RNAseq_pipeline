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
if (Org == 'hsa'){
    db <- 'org.Hs.eg.db'
    }
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

kegg_df <- enrichKEGG(target_gene_id_changed,
                      organism=Org,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      use_internal_data = T
                      )
kk = setReadable(kegg_df,
                  OrgDb = db,
                  keyType = "ENTREZID")
write.csv(data.frame(kk@result),file=outfile,quot=FALSE)

## draw
display_num <- 20
num <- nrow(as.data.frame(kegg_df))
if (num < display_num){
    #pdf(file=paste0(out_path,'/',prefix,'_kegg_enrichment_barplot.pdf'))
    #print(barplot(kegg_df, showCategory = num))
    #dev.off()
    pdf(file=paste0(out_path,'/',prefix,'_kegg_enrichment.pdf'))
    print(dotplot(kegg_df, showCategory = num))
    dev.off()
    }else{
    #pdf(file=paste0(out_path,'/',prefix,'_kegg_enrichment_barplot.pdf'))
    #print(barplot(kegg_df, showCategory = display_num))
    #dev.off()
    pdf(file=paste0(out_path,'/',prefix,'_kegg_enrichment.pdf'))
    print(dotplot(kegg_df, showCategory = display_num))
    dev.off()
}
