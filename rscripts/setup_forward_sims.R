root.dir<-"~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)

measure.times <- seq(100,1000,100) ## times to sample the input fitness landscape
Nreps <- 5 ## how many times to run each learned landscape (cpp) 

batch.file.path <- "ABM/validation/sim_forward/run.bat"
base.config.dir <- "ABM/config/"
input.sim.dir <- "ABM/output/"
output.sim.dir <- "ABM/validation/sim_forward/"
cpp.cmd <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/bin/Debug/ABM.exe"
validation.config.dir <- paste0(output.sim.dir,"config/")
validation.landscape.dir <- paste0(output.sim.dir,"landscape/")

dir.create(validation.landscape.dir)
dir.create(validation.config.dir)

sim.names <- list.files(input.sim.dir)
sim.names <- sim.names[grepl("v00",sim.names)]
output.sim.dirs <- paste0(output.sim.dir,sim.names)
lapply(output.sim.dirs,dir.create)
input.sim.dirs <- paste0(input.sim.dir,sim.names)

output.sim.dirs <- unlist(lapply(1:length(input.sim.dirs), function(i){
  paste0(output.sim.dirs[[i]],"/",list.files(input.sim.dirs[[i]]),"/")
}))


input.sim.dirs <- unlist(lapply(input.sim.dirs, function(di){
  paste0(di,"/",list.files(di),"/")
}))

setup_sims <- function(i){
  di <- input.sim.dirs[i]
  do <- output.sim.dirs[i]
  dir.create(do)
  
  x <- proc_sim(di,times=measure.times)
  x_opt <- optimsimple(x)
  lscape_opt <- optim_tps(x_opt)
  
  ## the following few lines are annoying scripty stuff to deal with setting up 
  ## the simulation filepaths etc. IDK how to simplify
  simname <- unlist(strsplit(di,"/"))
  # We are expecting to get something like this:
  #> simname
  #[1] "ABM"         "output"      "v00_rep_001" "00000"
  sim.id <- simname[grepl("v00",simname)]
  sim.rep <- tail(simname,1)
  
  base.config.path <- paste0(base.config.dir,sim.id,".txt")
  fitness.landscape.path <- paste0(getwd(),"/",validation.landscape.dir,sim.id,"_",sim.rep,".txt")
  config.path <- paste0(getwd(),"/",validation.config.dir,sim.id,"_",sim.rep,".txt")
  cpp.out.path <- paste0(getwd(),"/",do)
  config <- readLines(base.config.path)
  c.args <- data.frame(names=c("fitness_landscape_type","fitness_landscape_file","output_dir"),
                       vals=c("polyh",fitness.landscape.path,cpp.out.path))
  config <- modify_config(config,c.args)
  writeLines(config,config.path)
  fscape <- cbind(lscape_opt$knots,lscape_opt$w)
  fscape <- rbind(fscape,lscape_opt$v)
  write.table(fscape,fitness.landscape.path,row.names = FALSE,col.names=FALSE,sep=",")
  
  cmd <- paste(cpp.cmd,config.path)
  return(rep(cmd,Nreps))
}

cmds <- unlist(pbapply::pblapply(1:length(input.sim.dirs),setup_sims))
writeLines(cmds,batch.file.path)




