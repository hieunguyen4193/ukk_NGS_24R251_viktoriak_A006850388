gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"
library(clusterProfiler)
library(org.Hs.eg.db)
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

path.to.main.src.dir <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"

source(file.path(path.to.main.src.dir, "preparation.R"))
source(file.path(path.to.main.src.dir, "helper_functions_v2.R"))

condition.col <- "comparison2_WT"
condition1 <- "group1"
condition2 <- "group2"

outdir <- "/media/hieunguyen/HNSD01/outdir"

path.to.main.output <- file.path(outdir, "NGS_24R251_viktoriak_A006850388", "DESEQ2_output", condition.col)
path.to.save.output <- file.path(path.to.main.output, sprintf("case_%s_vs_%s", condition1, condition2))
path.to.save.pathway.output <- file.path(path.to.main.output, "pathway_analysis", sprintf("case_%s_vs_%s", condition1, condition2))
dir.create(path.to.save.pathway.output, showWarnings = FALSE, recursive = TRUE)

deseq.output <-  readRDS(file.path(path.to.save.output, "deseq_output.rds"))
resdf <- deseq.output$all.resultdf

convertdf <- bitr(resdf$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
resdf <- merge(resdf, convertdf, by.x = "gene_name", by.y = "SYMBOL")

#####--------------------------------------------------------------------#####
##### ORA WITH GO
#####--------------------------------------------------------------------#####
ora.GO <- enrichGO(gene = subset(resdf, resdf$padj <= 0.05)$gene_name,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE,
                   keyType = "SYMBOL",
                   pAdjustMethod = "BH")

if (is.null(ora.GO) == TRUE){
  ora.GOdf <- data.frame(status = c("No results from ORA - GO"))  
} else {
  ora.GOdf <- as.data.frame(ora.GO)
  ora.GOdf <- ora.GOdf[order(ora.GOdf$p.adjust, decreasing = FALSE),]  
}

dotplot(ora.GO, size = "GeneRatio", showCategory=30) + ggtitle("Dotplot for Over-represnetation analysis")

#####--------------------------------------------------------------------#####
##### ORA WITH KEGG
#####--------------------------------------------------------------------#####
ora.KEGG <-  enrichKEGG(gene = subset(resdf, resdf$padj <= 0.05)$ENTREZID,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05)
if (is.null(ora.KEGG) == TRUE){
  ora.KEGGdf <- data.frame(status = c("No results obtained from KEGG-ORA"))
} else {
  ora.KEGGdf <- as.data.frame(ora.KEGG)
  ora.KEGGdf <- ora.KEGGdf[order(ora.KEGGdf$p.adjust), ]    
}

#####--------------------------------------------------------------------#####
##### GSEA WITH GO
#####--------------------------------------------------------------------#####
tmp.full.list <- resdf %>% arrange(desc(log2FoldChange))
input.gene.list <- tmp.full.list$log2FoldChange
names(input.gene.list) <- tmp.full.list$gene_name

GSEA.GO <- gseGO(geneList = input.gene.list,
                 OrgD = org.Hs.eg.db,
                 ont = "ALL",
                 minGSSize = 100,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 keyType = "SYMBOL", seed = TRUE)

GSEA.GOdf <- as.data.frame(GSEA.GO) 

GSEA.GOdf <- GSEA.GOdf %>% rowwise() %>% 
  mutate(abs.NES = abs(NES)) %>%
  rownames_to_column("idx")

GSEA.GOdf <- GSEA.GOdf[order(GSEA.GOdf$NES, decreasing = TRUE), ]

top10_up_nes_pw.GO <- head(GSEA.GOdf, 10)$Description
top10_down_nes_pw.GO <- tail(GSEA.GOdf, 10)$Description

top10_up_nes_idx.GO <- head(GSEA.GOdf, 10)$idx
top10_down_nes_idx.GO <- tail(GSEA.GOdf, 10)$idx

#####--------------------------------------------------------------------#####
##### GSEA WITH KEGG
#####--------------------------------------------------------------------#####
tmp.full.list <- resdf %>% arrange(desc(log2FoldChange))
input.gene.list <- tmp.full.list$log2FoldChange
names(input.gene.list) <- convertdf$ENTREZID

GSEA.KEGG <- gseKEGG(geneList = input.gene.list,
                     organism = "hsa",
                     minGSSize = 100,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE, seed = TRUE)

GSEA.KEGGdf <- as.data.frame(GSEA.KEGG)

GSEA.KEGGdf <- GSEA.KEGGdf %>% rowwise() %>% 
  mutate(abs.NES = abs(NES)) %>%
  rownames_to_column("idx")

GSEA.KEGGdf <- GSEA.KEGGdf[order(GSEA.KEGGdf$NES, decreasing = TRUE), ]

top10_up_nes_pw.KEGG <- head(GSEA.KEGGdf, 10)$Description
top10_down_nes_pw.KEGG <- tail(GSEA.KEGGdf, 10)$Description

top10_up_nes_idx.KEGG <- head(GSEA.KEGGdf, 10)$idx
top10_down_nes_idx.KEGG <- tail(GSEA.KEGGdf, 10)$idx