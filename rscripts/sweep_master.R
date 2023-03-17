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
stop()
script_dir <- paste0(root.dir,"rscripts/")
cpp_source <- paste0(root.dir,"ABM/bin/Debug/ABM.exe")
if(root.dir)
sweep_dir <- paste0(root.dir,"ABM/output/sweep_01")
dir.create(sweep_dir)

source(paste0(script_dir,"sim_setup_functions.R"))
library(parallel)

fit_lscape <- function(setup_info){
  new.config.path <- tryCatch({
    source(setup_info$analfuncspath)
    dir <- setup_info$sweep_dir
    ntp <- setup_info$ntp
    ##perform the fitting
    sim_path <- paste0(dir,"train/00000/")
    measure.times=round(seq(3000/ntp,3000,3000/ntp))
    x <- proc_sim(sim_path,times=measure.times)
    x_opt <- optimsimple(x)
    x_opt$f_tru <- retrieve_fitness(rownames(x_opt), colnames(x$x), sim_path)
    optsimple.path <- paste0(dir,"optsimple_ntp_",ntp,".Rds")
    saveRDS(x_opt,optsimple.path)
    
    lscape_opt <- optim_tps(x_opt)
    
    ##make a config file (using original as basis) for the validation sim:
    base.config.path <- paste0(dir,"config.txt")
    new.config.path <- paste0(dir,"config_test.txt")
    cpp.out.path <- paste0(dir,"test")
    config <- readLines(base.config.path)
    new.landscape.path <- paste0(dir,"landscape_ntp_",ntp,".txt")
    c.args <- data.frame(names=c("fitness_landscape_type","fitness_landscape_file","output_dir"),
                         vals=c("polyh",new.landscape.path,cpp.out.path))
    config <- modify_config(config,c.args)
    writeLines(config,new.config.path)
    
    ##extract the fitted landscape and save it:
    fscape <- cbind(lscape_opt$knots,lscape_opt$w)
    fscape <- rbind(fscape,lscape_opt$v)
    write.table(fscape,new.landscape.path,row.names = FALSE,col.names=FALSE,sep=",")
    return(new.config.path)
  },error = function(e) return(NULL))
  
  
  return(new.config.path)
}

wavelength_range <- c(0.5,1,2,4)
Nchrom_range <- 2:12
ntp_range <- seq(4,20,2)


N_ls_reps <- 10
Ncores <- 3

cl <- makeCluster(getOption("cl.cores", Ncores))

for(Nchrom in Nchrom_range){
  for(wavelength in wavelength_range){
    setup_info <- lapply(1:N_ls_reps, function(dummyvar){
      gen_replicate(Nchrom=Nchrom,wavelength = wavelength,sweep_dir = sweep_dir,cpp_source = cpp_source)
    })
    
    parSapplyLB(cl,setup_info,function(xx) system(xx$cpp_run_cmd))
    setwd(sweep_dir)
    
    ##extra info that needs to get into the parallel function directly below:
    for(ntp in ntp_range){
      setup_info_ntp <- lapply(setup_info, function(xi){
        xi$sweep_dir <- paste0(sweep_dir,"/",xi$sweep_id,"/")
        xi$analfuncspath <- paste0(script_dir,"analysis_functions.R")
        xi$ntp=ntp
        return(xi)
      })
      
      new_config_paths <- parSapplyLB(cl,setup_info_ntp, fit_lscape)
      new_config_paths<-unlist(new_config_paths)
      if(sum(!is.null(new_config_paths))>0){
        cmds <- paste(cpp_source,new_config_paths)
        parSapplyLB(cl=cl,cmds,system)        
      }

    }
  }
}











