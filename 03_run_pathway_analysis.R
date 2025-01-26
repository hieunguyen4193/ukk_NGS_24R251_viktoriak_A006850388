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

outdir <- "/media/hieunguyen/HNSD01/outdir"
path.to.main.output <- file.path(outdir, "NGS_24R251_viktoriak_A006850388", "DESEQ2_output")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")

all.data <- Sys.glob(file.path(path.to.02.output, "*"))

all.cases <- basename(all.data)

for (condition.col in all.cases){
  all.comparisons <- Sys.glob(file.path(path.to.02.output, condition.col, "*"))
  for (input.compa in all.comparisons){
    condition1 <- str_split(basename(input.compa), "_")[[1]][[2]]
    condition2 <- str_split(basename(input.compa), "_")[[1]][[4]]
    print(sprintf("Working on condition %s, %s vs %s", condition.col, condition1, condition2))
    path.to.03.output <- file.path(path.to.main.output, "03_output", condition.col, sprintf("%s_vs_%s", condition1, condition2))
    dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
    
    deseq.output <-  readRDS(file.path(path.to.02.output, condition.col, sprintf("case_%s_vs_%s", condition1, condition2), "deseq_output.rds"))
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
    
    dotplot.ORA.GO <- dotplot(ora.GO, size = "GeneRatio", showCategory=30) + ggtitle("Dotplot for Over-represnetation analysis")
    writexl::write_xlsx(ora.GOdf, file.path(path.to.04.output, "ORA_GO_result.xlsx"))
    
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
    
    dotplot.ORA.KEGG <- dotplot(ora.KEGG, size = "GeneRatio", showCategory=30) + ggtitle("Dotplot for Over-represnetation analysis")
    writexl::write_xlsx(ora.KEGGdf, file.path(path.to.04.output, "ORA_KEGG_result.xlsx"))
    
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
    saveRDS(object = GSEA.GO, file.path(path.to.04.output, "GSEA.GO.rds"))
    writexl::write_xlsx(GSEA.GOdf, file.path(path.to.04.output, "GSEA_GOdf.xlsx"))
    
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
    
    saveRDS(object = GSEA.KEGG, file.path(path.to.04.output, "GSEA.KEGG.rds"))
    writexl::write_xlsx(GSEA.KEGGdf, file.path(path.to.04.output, "GSEA_KEGGdf.xlsx")) 
  }
}


