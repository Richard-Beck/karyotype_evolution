library(parallel)
library(fields)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"

script.dir <- paste0(root.dir,"rscripts/")
abm.out.dir <- paste0(root.dir,"ABM/output/sweep_01/")
sweep.dir <-"proc_data/02_testing_nn_inference/output/"
save.dir <- "proc_data/04_testing_spline_validation.Rds"

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
      
      cc <- sapply(1:3, function(dummyvar){
        xsam <- sample(1:nrow(xo),round(nrow(xo)*0.8))
        xtrain <- xmat[xsam,]
        ytrain <- y[xsam]
        xtest <- xmat[-xsam,]
        ytest<-y[-xsam]
        fit <- Krig(xtrain,ytrain,m=1,give.warnings = F)
        ypred <- predict(fit,xtest)
        cor=cor(ytest,ypred)
      })

      data.frame(nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],
                 r2=mean(cc*abs(cc)),model="optimsimple")
    },error=function(e) return(NULL))
  }
  if(which.model==3){
    a <- tryCatch({
      xo <- readRDS(paste0(root.dir,sweep.dir,"gfree/",f1,".Rds"))
      xmat <- do.call(rbind,lapply(rownames(xo), function(i){
        as.numeric(unlist(strsplit(i,split="[.]")))
      }))
      y <- xo$f_est
      
      cc <- sapply(1:3, function(dummyvar){
        xsam <- sample(1:nrow(xo),round(nrow(xo)*0.8))
        xtrain <- xmat[xsam,]
        ytrain <- y[xsam]
        xtest <- xmat[-xsam,]
        ytest<-y[-xsam]
        fit <- Krig(xtrain,ytrain,m=1,give.warnings = F)
        ypred <- predict(fit,xtest)
        cor=cor(ytest,ypred)
      })
      
      data.frame(nchrom=ids["N"],wl=ids["w"],rep=ids["rep"],
                 r2=mean(cc*abs(cc)),model="gfree")
    },error=function(e) return(NULL))
  }
  return(a)
}


setwd(abm.out.dir)
cl <- makeCluster(getOption("cl.cores", 2))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"analysis_functions_v2.R"))
  library(fields)
},script.dir=script.dir)
source(paste0(script.dir,"analysis_functions_v2.R"))
f0 <- list.files()#[1:5]
res <- list()
for(which.model in 2:3){
  x <- lapply(f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
  #x <- parLapplyLB(cl=cl,f0,optim_models,which.model=which.model,save.dir=save.dir,root.dir=root.dir,sweep.dir=sweep.dir)
  res <- c(res,list(do.call(rbind,x)))
}
res <- do.call(rbind,res)
saveRDS(res,paste0(root.dir,save.dir))











