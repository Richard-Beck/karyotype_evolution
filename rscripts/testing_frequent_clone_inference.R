library(parallel)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
script.dir <- paste0(root.dir,"rscripts/")
abm.out.dir <- paste0(root.dir,"ABM/output/nchrom_22/")
save.dir <- "proc_data/01_testing_frequent_clone_inference/"
dir.create(save.dir)
save.dir <- paste0(save.dir,"output/")
dir.create(save.dir)
dir.create(paste0(save.dir,"gfree"))
dir.create(paste0(save.dir,"optimsimple"))

optim_models <- function(id,which.model,save.dir,root.dir){
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
  xi <- proc_sim(dir,times=seq(400,2800,400))

  if(which.model==1){
    tryCatch({
      x <- optimise_frequent_clones(xi,min_obs = 20,use_gdata = T)
      pks <- lapply(rownames(x), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
      x$f_tru <- sapply(pks,getf,fobj=fobj)
      saveRDS(x,file=paste0(root.dir,save.dir,"optimsimple/",f1,".Rds"))
    },error=function(e) return(0))
  }
  if(which.model==2){
    tryCatch({
      x <- optimise_frequent_clones(xi,min_obs = 20,use_gdata = F)
      pks <- lapply(rownames(x), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
      x$f_tru <- sapply(pks,getf,fobj=fobj)
      saveRDS(x,file=paste0(root.dir,save.dir,"gfree/",f1,".Rds"))
    },error=function(e) return(0))
  }
  
}


setwd(abm.out.dir)
cl <- makeCluster(getOption("cl.cores", 3))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"analysis_functions_v2.R"))
},script.dir=script.dir)

f0 <- list.files()

for(which.model in 1:2){
  #sapply(f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir)
  parLapplyLB(cl=cl,f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir)
}











