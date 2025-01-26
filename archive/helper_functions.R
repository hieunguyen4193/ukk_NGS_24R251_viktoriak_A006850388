#######################################################################
# H E L P E R - F U N C T I O N S - DESEQ2 A N A L Y S I S 
#######################################################################

# trnguyen@ukaachen.de

# ---------------------------------------------------------------------
# generate DESeq2 dataset object
# ---------------------------------------------------------------------
generate_DESeq2_dataset <- function(path.to.salmon.quants, 
                                    path.to.tx2gene,
                                    metadata, 
                                    file.ext = "quant.sf",
                                    ERCC.spike.in = FALSE){
  # fetch all quant.sf files
  files <- file.path(path.to.salmon.quants, metadata$sample, file.ext)
  names(files) <- metadata$sample
  
  # read in the tx2gene list
  tx2gene <- read_tsv(path.to.tx2gene, 
                      col_names = c("transcript_id", "gene_id", "gene_name"), 
                      show_col_types = FALSE)
  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
  output <- DESeqDataSetFromTximport( txi,
                                      colData = metadata,
                                      design = ~ condition) # make sure that group exists in metadata
  
  ##### ##### ERCC SPIKE IN NORMALIZATION ##### #####
  
  ################################################### 
  return(output)
}

# ---------------------------------------------------------------------
# Run DESeq2 and preprocessing
# ---------------------------------------------------------------------
run_DESeq2_and_preprocess <- function(  deseq.dataset, 
                                        tx2gene, 
                                        thresh.pval = 0.05,
                                        geneid.to.remove = NULL){
  
  # run DESeq2
  if (is.null(geneid.to.remove) == FALSE){
    deseq.dataset <- deseq.dataset[rownames(deseq.dataset) %in% geneid.to.remove == FALSE]
  }
  
  deseq.object <- DESeq(deseq.dataset)
  
  # normalized read counts
  norm.counts <- counts(deseq.object, normalized=TRUE)
  
  # preprocessing
  ## extract the Differential expression analysis results from DESeq2 object
  ## this object contains: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  deseq.result <- results(deseq.object)
  
  # generate a dataframe containing gene.id and gene.name
  genedf <- data.frame(gene_id = tx2gene$gene_id, 
                       gene_name = tx2gene$gene_name) 
  
  ## remove duplicated gene.name / gene.id
  genedf <- genedf[!duplicated(genedf), ]
  
  ## merge the normalized count matrix and the gene information altogether
  norm.counts <- merge(genedf, norm.counts, by.x= "gene_id", by.y= "row.names")
  
  ## merge the normalized count matrix (with gene information) with DESeq2 results
  resultdf <- merge(genedf, as.data.frame(deseq.result), by.x= "gene_id", by.y="row.names", all.x=F, all.y=T)
  
  resultdf <- merge(resultdf, norm.counts, by=c("gene_id", "gene_name"), all.x=T, all.y=F)
  
  ## remove rows with NA values
  resultdf <- resultdf[complete.cases(resultdf), ]
  
  ## marking significant and non-significant genes by the DESeq2 p-value
  resultdf$sig <- "Non-sig."
  
  resultdf$sig[resultdf$padj < thresh.pval] <- "Sig. genes"
  
  ## Marking ERCC genes
  sel_ERCC <- str_detect(resultdf$gene_id, "^ERCC-*gene")
  
  resultdf$sig[sel_ERCC] <- "Spike in"
  
  ## remove ERCC spike-in genes
  resultdf <- resultdf[!sel_ERCC,]
  
  ## extract the dataframe with significant genes only
  resultdf.sig <- resultdf[resultdf$padj < thresh.pval,]
  
  ## extract the dataframe with significant genes only
  resultdf.nonsig <- resultdf[resultdf$padj >= thresh.pval,]
  
  # collect all outputs
  output <- list(norm.count = norm.counts,
                 all.resultdf = resultdf,
                 resultdf.sig = resultdf.sig,
                 resultdf.nonsig = resultdf.nonsig,
                 deseq.object = deseq.object)
  
  # return output
  return(output)
}

# ---------------------------------------------------------------------
# PCA PLOT
# ---------------------------------------------------------------------
plot.PCA.w.DESeq2.output <- function(input.df, 
                                     metadata,
                                     condition1,
                                     condition2){
  
  pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)
  
  pcadf <- data.frame(pca.object$x)
  
  row.names(pcadf) <- metadata$sample
  
  pcadf <- merge(pcadf, metadata, by.x = "row.names", by.y = "sample")
  
  pca.plot <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names)) +
    geom_point(size = 4) +
    ggtitle(sprintf("PCA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
  
  output = list(pca.df = pcadf,
                pca.plot = pca.plot)
  return(output)
}

# ---------------------------------------------------------------------
# B A S I C - A N A L Y S I S
# ---------------------------------------------------------------------
plot.heatmap <- function(selected.genes.to.heatmap, metadata, note, show.tick.labels = FALSE){
  input.to.heatmap <- subset(selected.genes.to.heatmap, select = c("gene_id", "gene_name", metadata$sample))
  
  heatmap.values <- log10(input.to.heatmap[,3:(dim(input.to.heatmap)[2])] + 1)
  
  selected.genes.heatmap.plot <- heatmaply(heatmap.values, 
                                           main= sprintf("Heatmap, Sample: %s vs. %s (%s)", condition1, condition2, note),
                                           method = "plotly",labRow=input.to.heatmap$gene_name,
                                           xlab = "Samples", ylab = "Genes", width = 800, height = 600,
                                           showticklabels = c(TRUE, show.tick.labels), show_dendrogram = c(FALSE, TRUE),
                                           key.title = "log10 scale colormap",
                                           label_names = c("Gene", "Sample", "Expression"),
                                           k_col = 2)
  
  return(selected.genes.heatmap.plot)
}

basic.analysis <- function( deseq.output, 
                            metadata, 
                            condition1,
                            condition2,
                            cutoff.logFC = 1, 
                            cutoff.adjp = 0.05){
  
  #######################################
  # PCA
  
  input.df <- deseq.output$norm.count[metadata$sample]
  
  pca.output <- plot.PCA.w.DESeq2.output(input.df = input.df, 
                                         metadata = metadata,
                                         condition1 = condition1, 
                                         condition2 = condition2)
  
  pca.plotly <- ggplotly(pca.output$pca.plot,  tooltip = c("all"))
  
  #######################################
  # volcano plot
  
  input.df <- deseq.output$all.resultdf
  
  input.df <- input.df %>% mutate(abs.log2FoldChange = abs(log2FoldChange))
  
  # input.df$show_label <- NA
  
  # condition on which gene_name should be shown, 
  # - significant level: adjusted p value less than cutoff.adjp
  # - log2FC greater than the cutoff.logFC
  
  # input.df$show_label[(input.df$abs.log2FoldChange > cutoff.logFC) & (input.df$padj <= cutoff.adjp)] <- input.df$gene_name[(input.df$abs.log2FoldChange > cutoff.logFC) & (input.df$padj <= cutoff.adjp)]
  
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=log2FoldChange, y=-log10(padj), col=sig, label=gene_name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    #geom_text() +
    ggtitle(sprintf("Volcano plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12))
  
  volcano.plotly <- ggplotly(volcano.plot,  tooltip = c("x", "y", "label"))
  #######################################
  # MA plot
  
  ma.plot <- ggplot(data=input.df, 
                    aes(x=log2(baseMean), y=log2FoldChange, col=sig, label=gene_name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    #geom_text() +
    ggtitle(sprintf("MA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12))
  
  ma.plotly <- ggplotly(ma.plot,  tooltip = c("x", "y", "label"))
  
  #######################################
  # output table on selected differentially expressed genes based on cutoff.logFC.
  
  sig.genes.with.highlogFC <- subset(input.df, (input.df$sig == "Sig. genes") & (input.df$abs.log2FoldChange > cutoff.logFC))
  
  nonsig.genes.with.highlogFC <- subset(input.df, (input.df$sig != "Sig. genes") & (input.df$abs.log2FoldChange > cutoff.logFC))
  
  #######################################
  # heatmap
  plot.heatmap <- function(selected.genes.to.heatmap, metadata, note, show.tick.labels = FALSE){
    input.to.heatmap <- subset(selected.genes.to.heatmap, select = c("gene_id", "gene_name", metadata$sample))
    
    heatmap.values <- log10(input.to.heatmap[,3:(dim(input.to.heatmap)[2])] + 1)
    
    selected.genes.heatmap.plot <- heatmaply(heatmap.values, 
                                             main= sprintf("Heatmap, Sample: %s vs. %s (%s)", condition1, condition2, note),
                                             method = "plotly",labRow=input.to.heatmap$gene_name,
                                             xlab = "Samples", ylab = "Genes", width = 800, height = 600,
                                             showticklabels = c(TRUE, show.tick.labels), show_dendrogram = c(FALSE, TRUE),
                                             key.title = "log10 scale colormap",
                                             label_names = c("Gene", "Sample", "Expression"),
                                             k_col = 2)
    
    return(selected.genes.heatmap.plot)
  }
  if (cutoff.logFC != 0){
    if (nrow(sig.genes.with.highlogFC) >= 2){
      sig.gene.heatmap.plotly <- plot.heatmap(sig.genes.with.highlogFC, metadata, "Sig.genes")
    } else {
      sig.gene.heatmap.plotly <- NA
    }
    if (nrow(nonsig.genes.with.highlogFC) >= 2){
      nonsig.gene.heatmap.plotly <- plot.heatmap(nonsig.genes.with.highlogFC, metadata, "Non-Sig.genes")
    } else {
      nonsig.gene.heatmap.plotly <- NA
    }
    
  }
  else{
    sig.gene.heatmap.plotly <- plot.heatmap(sig.genes.with.highlogFC, metadata, "Sig.genes")
    nonsig.gene.heatmap.plotly <- NA
  }
  
  # sig.gene.heatmap.plotly <- plot.heatmap(sig.genes.with.highlogFC, metadata, "Sig.genes")
  # nonsig.gene.heatmap.plotly <- plot.heatmap(nonsig.genes.with.highlogFC, metadata, "Non-Sig.genes")
  #######################################
  
  # generate output
  output <- list(
    sig.genes.with.highlogFC = sig.genes.with.highlogFC,
    nonsig.genes.with.highlogFC = nonsig.genes.with.highlogFC,
    pca.plot = pca.plotly,
    volcano.plot = volcano.plotly,
    ma.plot = ma.plotly,
    sig.gene.heatmap = sig.gene.heatmap.plotly,
    nonsig.gene.heatmap = nonsig.gene.heatmap.plotly
  )
  
}

# ---------------------------------------------------------------------
# function to generate table with interactive buttons in rmarkdown html
# ---------------------------------------------------------------------

create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

# ---------------------------------------------------------------------
# TRIWISE PLOT
# https://zouter.github.io/triwise/vignette.html 
# ---------------------------------------------------------------------






