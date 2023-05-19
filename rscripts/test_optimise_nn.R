library(parallel)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
script.dir <- paste0(root.dir,"rscripts/")
abm.out.dir <- paste0(root.dir,"ABM/output/nchrom_22/")
sweep.dir <-"proc_data/01_testing_frequent_clone_inference/output/"
save.dir <- "proc_data/02_testing_nn_inference/"
dir.create(save.dir)
save.dir <- paste0(save.dir,"output/")
dir.create(save.dir)
dir.create(paste0(save.dir,"gfree"))
dir.create(paste0(save.dir,"optimsimple"))

optim_models <- function(id,which.model,save.dir,root.dir,sweep.dir){
  print(id)
  f1 <- id
  id <- unlist(strsplit(id,"_"))
  ids <- id[c(2,4,6)]
  names(ids) <- id[c(1,3,5)]
  
  cfig_path <- paste0(f1,"/config.txt")
  lscape_path <- paste0(f1,"/landscape.txt")
  fobj <- gen_fitness_object(cfig_path, lscape_path)
  
  foi <- list.files(paste0(f1,"/train"))[1]
  dir <- paste0(f1,"/train/",foi)
  x <- tryCatch({proc_sim(dir,times=seq(400,2800,400))},
                error=function(e) return(NaN))
  a <- 0
  if(which.model==1){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"optimsimple/",f1,".Rds"))
      x2 <- wrap_neighbor_fitness(x,xo,use_gdata = T)
      pks <- lapply(rownames(x2), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
      x2$f_tru <- sapply(pks,getf,fobj=fobj)
      xo$stage <- "frequent"
      x2$stage <- "neighbor"
      xo <- rbind(xo,x2)
      saveRDS(xo,file=paste0(root.dir,save.dir,"optimsimple/",f1,".Rds"))
    },error=function(e) return(e))
  }
  if(which.model==2){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"gfree/",f1,".Rds"))
      x2 <- wrap_neighbor_fitness(x,xo,use_gdata = F)
      pks <- lapply(rownames(x2), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
      x2$f_tru <- sapply(pks,getf,fobj=fobj)
      xo$stage <- "frequent"
      x2$stage <- "neighbor"
      xo <- rbind(xo,x2)
      saveRDS(xo,file=paste0(root.dir,save.dir,"gfree/",f1,".Rds"))
    },error=function(e) return(e))
  }
  return(a)
}


setwd(abm.out.dir)
cl <- makeCluster(getOption("cl.cores", 3))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"analysis_functions_v2.R"))
},script.dir=script.dir)
source(paste0(script.dir,"analysis_functions_v2.R"))
f0 <- list.files()

for(which.model in 1:2){
  #x <- sapply(f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
  x <- parLapplyLB(cl=cl,f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
}











