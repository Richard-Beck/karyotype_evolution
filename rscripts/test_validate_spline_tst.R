
library(parallel)
library(fields)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
script.dir <- paste0(root.dir,"rscripts/")
abm.out.dir <- paste0(root.dir,"ABM/output/nchrom_22/")
save.dir <- "proc_data/04_testing_spline_validation/"
sweep.dir <- "proc_data/01_testing_frequent_clone_inference/output/"
dir.create(save.dir)

optim_models <- function(id,which.model,save.dir,root.dir,sweep.dir){
  print(id)
  f1 <- id
  id <- unlist(strsplit(id,"_"))
  ids <- id[c(2,4,6)]
  names(ids) <- id[c(1,3,5)]
  
  foi <- list.files(paste0(f1,"/train"))[1]
  dir <- paste0(f1,"/train/",foi)
  x <- proc_sim(dir,times=seq(400,2800,400))
  a <- NULL
  if(which.model==2){
    a <-tryCatch({
      
      xo <- readRDS(paste0(root.dir,sweep.dir,"optimsimple/",f1,".Rds"))
      x2 <- readRDS(paste0(root.dir,"proc_data/02_testing_nn_inference/output/optimsimple/",f1,".Rds"))
      xmat <- do.call(rbind,lapply(rownames(x2), function(i){
        as.numeric(unlist(strsplit(i,split="[.]")))
      }))
      ypred <- sapply(rownames(xo), function(xi){
        x2i <- xmat[!rownames(x2)==xi,]
        y2i <- x2$f_est[!rownames(x2)==xi]
        fit <- Krig(x2i,y2i)
        ytest <- matrix(as.numeric(unlist(strsplit(xi,split="[.]"))),nrow=1)
        predict(fit,ytest)
      })
      
      df <- data.frame(ytru=xo$f_tru,yclone=xo$f_est,ykrig=ypred,
                       nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],method="optimsimple",
                       row.names = rownames(xo))
      saveRDS(df,file=paste0(root.dir,save.dir,"optimsimple/",f1,".Rds"))
    },error=function(e) return(NULL))
  }
  if(which.model==3){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"gfree/",f1,".Rds"))
      x2 <- readRDS(paste0(root.dir,"proc_data/02_testing_nn_inference/output/gfree/",f1,".Rds"))
      xmat <- do.call(rbind,lapply(rownames(x2), function(i){
        as.numeric(unlist(strsplit(i,split="[.]")))
      }))
      ypred <- sapply(rownames(xo), function(xi){
        x2i <- xmat[!rownames(x2)==xi,]
        y2i <- x2$f_est[!rownames(x2)==xi]
        fit <- Krig(x2i,y2i)
        ytest <- matrix(as.numeric(unlist(strsplit(xi,split="[.]"))),nrow=1)
        predict(fit,ytest)
      })
      
      df <- data.frame(ytru=xo$f_tru,yclone=xo$f_est,ykrig=ypred,
                       nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],method="gfree",
                       row.names = rownames(xo))
      saveRDS(df,file=paste0(root.dir,save.dir,"gfree/",f1,".Rds"))
    },error=function(e) return(NULL))
  }
  return(a)
}


setwd(abm.out.dir)
cl <- makeCluster(getOption("cl.cores", 2))
clusterCall(cl, function(script.dir){
  library(fields)
  source(paste0(script.dir,"analysis_functions_v2.R"))
},script.dir=script.dir)
source(paste0(script.dir,"analysis_functions_v2.R"))
f0 <- list.files()

df <- parLapplyLB(cl=cl,f0[c(1,50,100,110)],optim_models,which.model=3,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
df <- do.call(rbind,df)
saveRDS(df,paste0(root.dir,save.dir,"gfree.Rds"))

df <- parLapplyLB(cl=cl,f0,optim_models,which.model=2,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
df <- do.call(rbind,df)
saveRDS(df,paste0(root.dir,save.dir,"optimsimple.Rds"))

#df <- lapply(f0,optim_models,which.model=3,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)












