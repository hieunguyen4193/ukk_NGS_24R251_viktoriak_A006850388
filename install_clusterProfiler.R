#####----------------------------------------------------------------------#####
##### UPGRADE THE CLUSTERPROFILER PACKAGES
devtools::install_github("YuLab-SMU/yulab.utils", upgrade = "never")
if ("clusterProfiler" %in% installed.packages() == TRUE){
  if (packageVersion("clusterProfiler") != "4.13.4"){
    remove.packages("DOSE")
    remove.packages("GOSemSim")
    remove.packages("yulab.utils")
    remove.packages("clusterProfiler")
    remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
    remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
  }
} else {
  remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
  remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
}
