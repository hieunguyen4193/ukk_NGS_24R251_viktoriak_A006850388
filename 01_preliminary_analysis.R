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
meta.data <- meta.data[, c("Sample.ID", "Name", "Genotype", "Treatment", condition.col)]
colnames(meta.data) <- c("sample", "Name", "Genotype", "Treatment", "condition")

data.samplesheet <- read.csv(file.path(path.to.main.src.dir, "SampleSheet.csv"))
outdir <- "/media/hieunguyen/HNSD01/outdir"

path.to.main.output <- file.path(outdir, "UKK_VK_RNAseq_data", "DESEQ2_output", condition.col)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
