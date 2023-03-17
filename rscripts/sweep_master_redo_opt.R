root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
script_dir <- paste0(root.dir,"rscripts/")
cpp_source <- paste0(root.dir,"ABM/bin/Debug/ABM.exe")
sweep_dir <- paste0(root.dir,"ABM/output/sweep_01")

source(paste0(script_dir,"sim_setup_functions.R"))
library(parallel)

fit_lscape <- function(setup_info){
  new.config.path <- tryCatch({
    options(scipen = 999)
    source(setup_info$analfuncspath)
    dir <- setup_info$sweep_dir
    ntp <- setup_info$ntp
    cfig_path <- paste0(dir,"/config.txt")
    lscape_path <- paste0(dir,"/landscape.txt")
    fobj <- gen_fitness_object(cfig_path, lscape_path)
    ##perform the fitting
    sim_path <- paste0(dir,"train/00000/")
    measure.times=round(seq(3000/ntp,3000,3000/ntp))
    x <- proc_sim(sim_path,times=measure.times)
    x_opt <- optimsimple(x)
    if(nrow(x_opt)<5) x_opt <- infer_nn(x,x_opt,pm0=0.00005)
    pks <- lapply(rownames(x_opt), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
    x_opt$f_tru <- sapply(pks,getf,fobj=fobj)
    optsimple.path <- paste0(dir,"optsimple_ntp_",ntp,".Rds")
    saveRDS(x_opt,optsimple.path)
    
    lscape_opt <- optim_tps(x_opt)
    
    ##make a config file (using original as basis) for the validation sim:
    base.config.path <- paste0(dir,"config.txt")
    new.config.path <- paste0(dir,"config_ntp_",ntp,".txt")
    cpp.out.path <- paste0(dir,"test5")
    config <- readLines(base.config.path)
    new.landscape.path <- paste0(dir,"t5_landscape_ntp_",ntp,".txt")
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

setwd(sweep_dir)

ids <- list.files()

wavelength_range <- c(0.5,1,2,4)
Nchrom_range <- 2:12
ntp_range <- seq(4,12,2)


N_ls_reps <- 10
Ncores <- 3

cl <- makeCluster(getOption("cl.cores", Ncores))

for(ntp in ntp_range){
  setup_info_ntp <- lapply(ids, function(xi){
    dir.create(paste0(xi,"/test5"))
    list(cpp_run_cmd=cpp_source, sweep_id=xi,ntp=ntp,
         sweep_dir = paste0(sweep_dir,"/",xi,"/"),
         analfuncspath = paste0(script_dir,"analysis_functions.R"))
  })
  #new_config_paths <- sapply(setup_info_ntp, fit_lscape)
  new_config_paths <- parSapplyLB(cl,setup_info_ntp, fit_lscape)
  new_config_paths<-unlist(new_config_paths)
  if(sum(!is.null(new_config_paths))>0){
    cmds <- paste(cpp_source,new_config_paths)
    parSapplyLB(cl=cl,cmds,system)        
  }
}








