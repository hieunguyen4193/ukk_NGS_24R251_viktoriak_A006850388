gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

path.to.main.src.dir <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"

source(file.path(path.to.main.src.dir, "preparation.R"))
source(file.path(path.to.main.src.dir, "helper_functions_v2.R"))

path.to.metadata <- file.path(path.to.main.src.dir, "metadata.csv")
meta.data <- read.csv(path.to.metadata)
meta.data <- meta.data[, c("Sample.ID", "Name", "Genotype", "Treatment", "batch")]
colnames(meta.data) <- c("sample", "Name", "Genotype", "Treatment", "condition")

outdir <- "/media/hieunguyen/HNSD01/outdir"

path.to.main.output <- file.path(outdir, "NGS_24R251_viktoriak_A006850388", "DESEQ2_output")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.salmon.output <- file.path(outdir, "NGS_24R251_viktoriak_A006850388", "star_salmon")
path.to.tx2gene <- file.path(path.to.salmon.output, "tx2gene.tsv")

deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = path.to.salmon.output,
                                         path.to.tx2gene = path.to.tx2gene,
                                         metadata = meta.data)

norm.count.mat <- vst(deseq.dataset, blind = FALSE)

pcaData <- plotPCA(norm.count.mat, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca.plot.batch <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
