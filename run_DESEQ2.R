gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"
  
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

#####----------------------------------------------------------------------#####
##### DEFINE FUNCTION TO RUN DESEQ2
#####----------------------------------------------------------------------#####
run_pipeline_DESEQ2 <- function(condition1, condition2, condition.col){
  path.to.main.src.dir <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"
  
  source(file.path(path.to.main.src.dir, "preparation.R"))
  source(file.path(path.to.main.src.dir, "helper_functions_v2.R"))
  
  path.to.metadata <- file.path(path.to.main.src.dir, "metadata.csv")
  meta.data <- read.csv(path.to.metadata)
  meta.data <- meta.data[, c("Sample.ID", "Name", "Genotype", "Treatment", condition.col)]
  colnames(meta.data) <- c("sample", "Name", "Genotype", "Treatment", "condition")
  
  data.samplesheet <- read.csv(file.path(path.to.main.src.dir, "SampleSheet.csv"))
  outdir <- "/media/hieunguyen/HNSD01/outdir"
  
  path.to.main.output <- file.path(outdir, "NGS_24R251_viktoriak_A006850388", "DESEQ2_output", condition.col)
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.save.output <- file.path(path.to.main.output, sprintf("case_%s_vs_%s", condition1, condition2))
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  
  filtered.metadata <- subset(meta.data, meta.data$condition %in% c(condition1, condition2))
  filtered.metadata$condition <- factor(filtered.metadata$condition, levels = c(condition1, condition2))
  
  salmon.output <- "/media/hieunguyen/HNSD01/outdir/NGS_24R251_viktoriak_A006850388/star_salmon"
  
  ##### save metadata
  writexl::write_xlsx(filtered.metadata, file.path(path.to.save.output, "metadata_for_this_analysis.xlsx"))
  
  ##### DESEQ2
  path.to.tx2gene <- file.path(salmon.output, "tx2gene.tsv")
  deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                           path.to.tx2gene = path.to.tx2gene,
                                           metadata = filtered.metadata)
  tx2gene <- read_tsv(path.to.tx2gene, 
                      col_names = c("transcript_id", "gene_id", "gene_name"), 
                      show_col_types = FALSE)
  deseq.output <- run_DESeq2_and_preprocess(deseq.dataset, tx2gene, thresh.pval = 0.05)
  saveRDS(deseq.output, file.path(path.to.save.output, "deseq_output.rds"))
  
  sigdf <- deseq.output$resultdf.sig
  nonsigdf <- deseq.output$resultdf.nonsig
  
  writexl::write_xlsx(sigdf, file.path(path.to.save.output, "sig_DGE.xlsx"))
  writexl::write_xlsx(nonsigdf, file.path(path.to.save.output, "nonsig_DGE.xlsx"))
  
  ##### PCA plot
  input.df <- deseq.output$norm.count[filtered.metadata$sample]
  pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)
  
  pcadf <- data.frame(pca.object$x)
  
  row.names(pcadf) <- filtered.metadata$sample
  
  pcadf <- merge(pcadf, filtered.metadata, by.x = "row.names", by.y = "sample")
  pcadf <- pcadf %>% rowwise() %>%
    mutate(full.name = sprintf("%s_%s_%s_%s", Row.names, Name, Genotype, Treatment), )
  
  pca.plot <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
    geom_point(size = 4) +
    ggtitle(sprintf("PCA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    geom_label_repel() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
  
  ggsave(plot = pca.plot, filename = "PCA.svg", path = path.to.save.output, dpi = 300, width = 10, height = 10)
  
  ##### Volcano plot
  cutoff.adjp <- 0.05
  cutoff.logFC <- 1
  
  input.df <- deseq.output$all.resultdf
  input.df <- input.df %>% mutate(abs.log2FoldChange = abs(log2FoldChange))
  input.df <- input.df %>% rowwise() %>%
    mutate(show.gene.name = ifelse(padj < cutoff.adjp, gene_name, NA))
  
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=log2FoldChange, y=-log10(padj), col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    geom_text_repel() +
    ggtitle(sprintf("Volcano plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    xlim(c(-max(input.df$abs.log2FoldChange), max(input.df$abs.log2FoldChange)))
  ggsave(plot = volcano.plot, filename = "volcano_plot.svg", path = path.to.save.output, dpi = 300, width = 14, height = 10)
  
  ##### MA plot
  ma.plot <- ggplot(data=input.df, 
                    aes(x=log2(baseMean), y=log2FoldChange, col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_text_repel() +
    ggtitle(sprintf("MA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12))
  ggsave(plot = ma.plot, filename = "MA_plot.svg", path = path.to.save.output, dpi = 300, width = 14, height = 10)
  
  ##### heatmap
  sig.genes.with.highlogFC <- subset(input.df, (input.df$sig == "Sig. genes") & (input.df$abs.log2FoldChange > cutoff.logFC))
  nonsig.genes.with.highlogFC <- subset(input.df, (input.df$sig != "Sig. genes") & (input.df$abs.log2FoldChange > cutoff.logFC))
  input.to.heatmap <- subset(sig.genes.with.highlogFC, select = c("gene_id", "gene_name", filtered.metadata$sample))
  if (nrow(input.to.heatmap) > 0){
    heatmap.values <- log10(input.to.heatmap[,3:(dim(input.to.heatmap)[2])] + 1)
    selected.genes.heatmap.plot <- heatmaply(heatmap.values, 
                                             main= sprintf("Heatmap, Sample: %s vs. %s", condition1, condition2),
                                             method = "plotly",labRow=input.to.heatmap$gene_name,
                                             xlab = "Samples", ylab = "Genes", width = 800, height = 600,
                                             showticklabels = c(TRUE, FALSE), show_dendrogram = c(FALSE, TRUE),
                                             key.title = "log10 scale colormap",
                                             label_names = c("Gene", "Sample", "Expression"),
                                             k_col = 2, file = file.path(path.to.save.output, "heatmap.html"))
  }
  
  ##### TK genes
  tk.genes <- read.csv(file.path(path.to.main.src.dir, "geneSet/TK_genes.csv"))
  new.gene.sets <- read.csv(file.path(path.to.main.src.dir, "geneSet/new_gene_set_20240409.csv"))
  plot.genes <- c(tk.genes$TK_gene, new.gene.sets$Gene) %>% unique()
  # mouse.genes <- to_vec(
  #   for (item in plot.genes){
  #     paste0(str_split(item, "")[[1]][[1]], 
  #            paste0(str_split(item, "")[[1]][2:nchar(item)] %>% tolower(), collapse = ""))
  #   }
  # )
  countdf <- deseq.output$norm.count %>% subset(select = -c(gene_id)) %>%
    subset(gene_name %in% plot.genes)
  countdf <- countdf[, c("gene_name", filtered.metadata$sample)]
  # countdf <- countdf[!duplicated(countdf$gene_name),]
  heatmaply(countdf %>% subset(select = -c(gene_name)),
            scale = "column",
            main= sprintf("Heatmap, Sample: %s vs. %s", condition1, condition2),
            method = "plotly",labRow=countdf$gene_name,
            xlab = "Samples", ylab = "Genes", width = 1500, height = 2400,
            showticklabels = c(TRUE, TRUE), show_dendrogram = c(FALSE, TRUE),
            key.title = "log10 scale colormap",
            label_names = c("Gene", "Sample", "Expression"),
            k_col = 2, file = file.path(path.to.save.output, "heatmap_TK_genes.html"))
}

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####

condition1 <- "group1"
condition2 <- "group2"
all.conditions <- c("comparison1",
                    "comparison2_WT",
                    "comparison2_CSK_KO",
                    "comparison2_LYN_KO")
for (condition.col in all.conditions){
  run_pipeline_DESEQ2(condition1 = "group1", 
                      condition2 = "group2",
                      condition.col = condition.col)
}
run_pipeline_DESEQ2(condition1 = "group1", 
                    condition2 = "group3",
                    condition.col = "comparison1")
run_pipeline_DESEQ2(condition1 = "group2", 
                    condition2 = "group3",
                    condition.col = "comparison1")


