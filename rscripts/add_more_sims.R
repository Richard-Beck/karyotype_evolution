root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
script_dir <- paste0(root.dir,"rscripts/")
cpp_source <- paste0(root.dir,"ABM/bin/Debug/ABM.exe")
sweep_dir <- paste0(root.dir,"ABM/output/sweep_01/")


runfunc <- function(id,root.dir,Nreps=5){
  cpp_source <- paste0(root.dir,"ABM/bin/Debug/ABM.exe")
  sweep_dir <- paste0(root.dir,"ABM/output/sweep_01/")
  cfig <- paste0(sweep_dir,id,"/config.txt")
  cmd <- paste(cpp_source,cfig)
  for(i in 1:Nreps){
    tryCatch(system(cmd),error=function(e) return(0))  
  }
}

library(parallel)
ids <- list.files(sweep_dir)
cl <- makeCluster(getOption("cl.cores", 3))
parLapplyLB(cl,ids,runfunc,root.dir=root.dir)


