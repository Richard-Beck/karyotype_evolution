source("~/projects/008_birthrateLandscape/karyotype_evolution/rscripts/landscape_functions.R")
source("~/projects/008_birthrateLandscape/karyotype_evolution/rscripts/analysis_functions_v2.R")

setwd("~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/sweep_01/")
f0 <- list.files()


lapply(f0,function(id){
  print(id)
  f1 <- id
  id <- unlist(strsplit(id,"_"))
  ids <- id[c(2,4,6)]
  names(ids) <- id[c(1,3,5)]
  
  cfig_path <- paste0(f1,"/config.txt")
  lscape_path <- paste0(f1,"/landscape.txt")
  fobj <- gen_fitness_object(cfig_path, lscape_path)
  
  fo <- list.files(paste0(f1,"/train"))[1:2]
  lapply(fo, function(foi){
    dir <- paste0(f1,"/train/",foi)
    xi <- proc_sim(dir,times=seq(400,2800,400))
    tryCatch({
      stop()
      x_gfree <- optimise_frequent_clones(xi,min_obs = 5,use_gdata = F)
      pks <- lapply(rownames(x_gfree), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
      x_gfree$f_tru <- sapply(pks,getf,fobj=fobj)
      saveRDS(x_gfree,file=paste0(f1,"/xopt_gfree_train_",foi,".Rds"))
    },error=function(e) return(0))

    tryCatch({
      x_gdata <- optimise_frequent_clones(xi,min_obs = 5,use_gdata = T)
      pks <- lapply(rownames(x_gdata), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
      x_gdata$f_tru <- sapply(pks,getf,fobj=fobj)
      saveRDS(x_gdata,file=paste0(f1,"/xopt_gdata3_train_",foi,".Rds"))
    },error=function(e) return(0))
    

    
  })
})










