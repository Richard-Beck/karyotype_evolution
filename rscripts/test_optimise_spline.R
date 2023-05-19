library(parallel)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
input <- "stage_1"
script.dir <- paste0(root.dir,"rscripts/")
abm.out.dir <- paste0(root.dir,"ABM/output/nchrom_22/")
sweep.dir <-"proc_data/02_testing_nn_inference/output/"
save.dir <- "proc_data/03_testing_spline_inference_nnin.Rds"

if(input=="stage_1"){
  sweep.dir <-"proc_data/01_testing_frequent_clone_inference/output/"
  save.dir <- "proc_data/03_testing_spline_inference_fqin.Rds"
}

optim_models <- function(id,which.model,save.dir,root.dir,sweep.dir){
  print(id)
  f1 <- id
  id <- unlist(strsplit(id,"_"))
  ids <- id[c(2,4,6)]
  names(ids) <- id[c(1,3,5)]
  
  cfig_path <- paste0(f1,"/config.txt")
  lscape_path <- paste0(f1,"/landscape.txt")
  fobj <- gen_fitness_object(cfig_path, lscape_path)
  
  a <- 0
  if(which.model==1){
    tryCatch({
      stop() ## previously tested another inference method but it was bad so abandoned it.
    },error=function(e) return(0))
  }
  if(which.model==2){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"optimsimple/",f1,".Rds"))
      xmat <- do.call(rbind,lapply(rownames(xo), function(i){
        as.numeric(unlist(strsplit(i,split="[.]")))
      }))
      y <- xo$f_est
      fit <- Krig(xmat,y,m=1,give.warnings = F)
      tst <- gen_all_neighbours(rownames(xo))
      if(!"stage"%in%colnames(xo)){
        tst <- lapply(1:nrow(tst),function(i) tst[i,])
        tst <- gen_all_neighbours(tst,as.strings=F)
        tst_names <- apply(tst,1,paste,collapse=".")
        tst <- tst[!tst_names%in%rownames(xo),]
      }
      #tst <- expand_pop(rownames(xo),fc=2000/nrow(xo))
      ytru <- apply(tst,1,getf,fobj=fobj)
      ypred <- predict(fit,tst)
      ytru<-ytru-mean(ytru)
      ypred<-ypred-mean(ypred)
      rsq <- R2(obs=ytru,pred=ypred)
      data.frame(nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],
                 cor=cor(ytru,ypred),R2=rsq,model="optimsimple",ntest=nrow(tst),n0=nrow(xo))
    },error=function(e) return(NULL))
  }
  if(which.model==3){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"gfree/",f1,".Rds"))
      xmat <- do.call(rbind,lapply(rownames(xo), function(i){
        as.numeric(unlist(strsplit(i,split="[.]")))
      }))
      y <- xo$f_est
      fit <- Krig(xmat,y,m=1,give.warnings = F)
      tst <- gen_all_neighbours(rownames(xo))
      if(!"stage"%in%colnames(xo)){
        tst <- lapply(1:nrow(tst),function(i) tst[i,])
        tst <- gen_all_neighbours(tst,as.strings=F)
        tst_names <- apply(tst,1,paste,collapse=".")
        tst <- tst[!tst_names%in%rownames(xo),]
      }
      
      #tst <- expand_pop(rownames(xo),fc=2000/nrow(xo))
      ytru <- apply(tst,1,getf,fobj=fobj)
      ypred <- predict(fit,tst)
      ytru<-ytru-mean(ytru)
      ypred<-ypred-mean(ypred)
      rsq <- R2(obs=ytru,pred=ypred)
      data.frame(nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],
                 cor=cor(ytru,ypred),R2=rsq,model="gfree",ntest=nrow(tst),n0=nrow(xo))
    },error=function(e) return(NULL))
  }
  print(a)
  return(a)
}


setwd(abm.out.dir)
cl <- makeCluster(getOption("cl.cores", 2))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"analysis_functions_v2.R"))
  library(fields)
},script.dir=script.dir)
source(paste0(script.dir,"analysis_functions_v2.R"))
f0 <- list.files()#[1:3]#[1:5]

res <- list()
for(which.model in 2:3){
  x <- lapply(f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
  #x <- parLapplyLB(cl=cl,f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
  res <- c(res,list(do.call(rbind,x)))
}
res <- do.call(rbind,res)
print(head(res))
saveRDS(res,paste0(root.dir,save.dir))











