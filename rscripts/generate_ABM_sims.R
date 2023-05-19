## https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

root.dir <- gsub("[\\]","/",thisFile())
root.dir <- unlist(strsplit(root.dir,split="/"))
root.dir <- root.dir[1:(length(root.dir)-2)]
root.dir <- paste0(paste(root.dir,collapse="/"),"/")

script_dir <- paste0(root.dir,"rscripts/")
cpp_source <- paste0(root.dir,"ABM/bin/Debug/ABM")
if(Sys.info()["sysname"]=="Windows") cpp_source <- paste0(cpp_source,".exe")

print("enter name for new sweep:")
sweep_name = readLines(con = "stdin", n = 1)

print("enter number of cores:")
Ncores = as.numeric(readLines(con = "stdin", n = 1))

sweep_dir <- paste0(root.dir,"ABM/output/",sweep_name)
dir.create(sweep_dir)

source(paste0(script_dir,"sim_setup_functions.R"))
library(parallel)

N_ls_reps <- 25
wavelength_range <- c(0.1,0.2,0.4,0.8,1.6)
cl <- makeCluster(getOption("cl.cores", min(N_ls_reps,Ncores)))
Nchrom <- 22
  for(wavelength in wavelength_range){
    setup_info <- lapply(1:N_ls_reps, function(dummyvar){
      gen_replicate(Nchrom=Nchrom,wavelength = wavelength,sweep_dir = sweep_dir,cpp_source = cpp_source)
    })
    parSapplyLB(cl,setup_info,function(xx) system(xx$cpp_run_cmd))
    
}












