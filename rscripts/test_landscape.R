cellLine <- "SA609"
base.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(base.dir)
x_opt <- readRDS(paste0("data/",cellLine,"/optim_stage1.Rds"))
options(scipen = 999)
source("rscripts/analysis_functions.R")
root.dir <- paste0(base.dir,"data/",cellLine,"/")
base_config_path <- paste0(root.dir,"ABM_inputs/krig_config.txt")
fitness.landscape.path <- paste0(root.dir,"ABM_inputs/krig_lscape.txt")
Nreps <- 3
Nstart <- 500000
base_out <- paste0(root.dir,"ABM_output/")
input.times <- 6:9

output_dirs <- paste0(base_out,"tstart_",input.times,"/")

sapply(output_dirs,dir.create)

base_config <- readLines(base_config_path)


c.args <- data.frame(names=c("bf","dt","p","init_size","max_size","Nsteps","fitness_landscape_file"),
                     vals=c(0.5,0.1,0.00005,Nstart,1000000,1000,fitness.landscape.path))
base_config <- modify_config(base_config,c.args)

x <- readRDS(paste0("data/",cellLine,"/opt_data.Rds"))

for(i in 1:length(input.times)){
  Tstart <- input.times[i]
  xi <- x$x[x$x[,Tstart]>3,Tstart]
  #if(input.times[i]==4){
   # xi <- x_opt$x0
    #names(xi) <- rownames(x_opt)
  #}
  xi <- round(xi*Nstart/sum(xi))
  xi <- xi[xi>0]
  xmat <- do.call(rbind,lapply(names(xi), function(ri){
    as.numeric(unlist(strsplit(ri,split="[.]")))
  }))
  xmat <- cbind(xmat,as.numeric(xi))
  xmat <- xmat[apply(xmat,1,function(xj) sum(xj==0)==0),]
  
  pop_filepath <- paste0(output_dirs[i],"popfile.txt")
  
  write.table(xmat,pop_filepath,row.names = FALSE,col.names=FALSE,sep=",")
  
  
  ABM_out_dir <- paste0(output_dirs[i],"ABM_output")
  dir.create(ABM_out_dir)
  
  cfigi <- base_config  
  c.args <- data.frame(names=c("population_file","output_dir"),
                       vals=c(pop_filepath,ABM_out_dir))
  cfigi <- modify_config(cfigi,c.args)
  
  cfig_path <- paste0(output_dirs[i],"config.txt")
  writeLines(cfigi,cfig_path)
  
  cmd <- paste("C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/bin/Debug/ABM.exe",cfig_path)
  for(i in 1:Nreps) system(cmd)
}






