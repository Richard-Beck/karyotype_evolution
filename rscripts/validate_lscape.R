library(fields)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
script.dir <- paste0(root.dir,"rscripts/")
source(paste0(script.dir,"analysis_functions_v2.R"))
setwd(root.dir)
optim_loo <- function(i,xx,xo){
  xi <- rownames(xo)[i]
  xistr <- xi
  print(xi)
  xi <-as.numeric(unlist(strsplit(xi,split="[.]"))) 
  xtrial <- xo
  xtrial <- xtrial[!rownames(xtrial)==xistr,]
  x2 <- wrap_neighbor_fitness(xx,xtrial,use_gdata = F)
  x2 <- x2[!rownames(x2)==xistr,]
  x2 <- rbind(xtrial,x2)
  xmat <- do.call(rbind,lapply(rownames(x2), function(i){
    as.numeric(unlist(strsplit(i,split="[.]")))
  }))
  y <- x2$f_est
  fit <- Krig(xmat,y,m=1)
  pred <- predict(fit,matrix(xi,nrow=1))
  data.frame(pred=pred,row.names = xistr)
}

validate_lscape <- function(cellLine){
  xx <- readRDS(paste0("salehi_data/02_optim_data/",cellLine,".Rds"))
  xx$dt <- 25
  if(grepl("hTERT",cellLine)) xx$dt <- 5
  xo <- readRDS(paste0("salehi_data/03_inference_output/",cellLine,"/frequent_clones.Rds"))
  df <- do.call(rbind,lapply(1:nrow(xo),function(i) optim_loo(i,xx,xo)))
  xo$validation <- df$pred
  saveRDS(xo,paste0("salehi_data/03_inference_output/",cellLine,"/validation.Rds"))
  print(R2(xo$f_est,xo$validation))
}

ff <- list.files("salehi_data/03_inference_output/")
cellLines <- sapply(ff,function(fi) head(unlist(strsplit(fi,split=".Rds")),1))
lapply(cellLines,validate_lscape)

