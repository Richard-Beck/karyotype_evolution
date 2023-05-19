
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
  optim_loo <- function(xi,xx,xo,use_gdata){
    xistr <- xi
    print(xi)
    xi <-as.numeric(unlist(strsplit(xi,split="[.]"))) 
    xx$x <- xx$x[,order(as.numeric(colnames(xx$x)))]
    #xx$x <- xx$x[-i,]
    
    #xo <- optimise_frequent_clones(xx,min_obs=20,use_gdata = use_gdata)
    xo <- xo[!rownames(xo)==xistr,]
    x2 <- wrap_neighbor_fitness(xx,xo,use_gdata = use_gdata)
    x2 <- x2[!rownames(x2)==xistr,]
    xo <- rbind(xo[,colnames(x2)],x2)
    xmat <- do.call(rbind,lapply(rownames(xo), function(i){
      as.numeric(unlist(strsplit(i,split="[.]")))
    }))
    y <- xo$f_est
    fit <- Krig(xmat,y,m=1)
    pred <- predict(fit,matrix(xi,nrow=1))
    list(xo=xo,xi=xi,xistr=xistr,pred=pred)
  }
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
      xtst <- lapply(rownames(xo), optim_loo,xx=x,xo=xo,use_gdata=T)
      ypred <- sapply(xtst,function(xi) xi$pred)
      
      df <- data.frame(ytru=xo$f_tru,yclone=xo$f_est,ykrig=ypred,
                 nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],method="optimsimple",
                 row.names = rownames(xo))
      saveRDS(df,file=paste0(root.dir,save.dir,"optimsimple/",f1,".Rds"))
    },error=function(e) return(NULL))
  }
  if(which.model==3){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"gfree/",f1,".Rds"))
      xtst <- lapply(rownames(xo), optim_loo,xx=x,xo=xo,use_gdata=F)
      ypred <- sapply(xtst,function(xi) xi$pred)
      
      df <- data.frame(ytru=xo$f_tru,yclone=xo$f_est,ykrig=ypred,
                       nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],method="gfree",
                       row.names = rownames(xo))
      saveRDS(df,file=paste0(root.dir,save.dir,"gfree/",f1,".Rds"))
    },error=function(e) return(NULL))
  }
  return(a)
}


setwd(abm.out.dir)
cl <- makeCluster(getOption("cl.cores", 50))
clusterCall(cl, function(script.dir){
  library(fields)
  source(paste0(script.dir,"analysis_functions_v2.R"))
},script.dir=script.dir)
source(paste0(script.dir,"analysis_functions_v2.R"))
f0 <- list.files()

#df <- parLapplyLB(cl=cl,f0,optim_models,which.model=2,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
#df <- do.call(rbind,df)
#saveRDS(df,paste0(root.dir,save.dir,"optimsimple.Rds"))

#df <- lapply(f0,optim_models,which.model=3,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
df <- parLapplyLB(cl=cl,f0,optim_models,which.model=3,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
df <- do.call(rbind,df)
saveRDS(df,paste0(root.dir,save.dir,"gfree.Rds"))









