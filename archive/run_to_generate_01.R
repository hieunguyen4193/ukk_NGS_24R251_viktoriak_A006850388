gc()
rm(list = ls())

path.to.main.src.dir <- "/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388"

source(file.path(path.to.main.src.dir, "preparation.R"))
source(file.path(path.to.main.src.dir, "helper_functions.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"
path.to.main.output <- file.path(outdir, "NGS_24R251_viktoriak_A006850388", "DESEQ2_output")
path.to.save.html <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html,showWarnings = FALSE, recursive = TRUE)

path.to.metadata <- file.path(path.to.main.src.dir, "metadata.csv")
meta.data <- read.csv(path.to.metadata)

path.to.rmd <- file.path(path.to.main.src.dir, "01_RNAseq_data_analysis.Rmd")

all.conditions <- c("comparison1",
                    "comparison2_WT",
                    "comparison2_CSK_KO",
                    "comparison2_LYN_KO")

for (condition.col in all.conditions){
  save.html.name <- sprintf("DESEQ2_analysis_%s.html", condition.col)
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(
      input = path.to.rmd, 
      params = list(
        condition1 = "group1",
        condition2 = "group2",
        condition.col = condition.col
      ),
      output_file = save.html.name,
      output_dir =  path.to.save.html
    )
  }
}
