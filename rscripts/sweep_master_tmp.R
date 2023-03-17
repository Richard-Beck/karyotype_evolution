root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
script_dir <- paste0(root.dir,"rscripts/")
cpp_source <- paste0(root.dir,"ABM/bin/Debug/ABM.exe")
sweep_dir <- paste0(root.dir,"ABM/output/sweep_01")

source(paste0(script_dir,"sim_setup_functions.R"))
source(paste0(script_dir,"analysis_functions.R"))

dirs <- list.files(sweep_dir)
ntp_range <- seq(4,20,2)

extract_optsimple <- function(dir){
  for(ntp in ntp_range){
    tryCatch({sim_path <- paste0(sweep_dir,"/",dir,"/train/00000/")
    measure.times=round(seq(3000/ntp,3000,3000/ntp))
    x <- proc_sim(sim_path,times=measure.times)
    x_opt <- optimsimple(x)
    x_opt$f_tru <- retrieve_fitness(rownames(x_opt), colnames(x$x), sim_path)
    optsimple.path <- paste0(sweep_dir,"/",dir,"/optsimple_ntp_",ntp,".Rds")
    saveRDS(x_opt,optsimple.path)},error=function(e) return(0))
    
  }
  return(0)
}

pbapply::pblapply(dirs,extract_optsimple)









