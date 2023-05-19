library(fields)
root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
script.dir <- paste0(root.dir,"rscripts/")
source(paste0(script.dir,"analysis_functions_v2.R"))
setwd(root.dir)
landscape_inference <- function(cellLine){
  print(cellLine)
  x <- readRDS(paste0("salehi_data/02_optim_data/",cellLine,".Rds"))
  x$x <- x$x[,order(as.numeric(colnames(x$x)))]
  x$dt <- 25
  if(grepl("hTERT",cellLine)) x$dt <- 5
  #diploid <- paste(rep(2,22),collapse=".")
  #x$x <- x$x[!rownames(x$x)==diploid,]
  xo <- optimise_frequent_clones(x,min_obs=10)
  dir.create(paste0("salehi_data/03_inference_output/",cellLine))
  saveRDS(xo,paste0("salehi_data/03_inference_output/",cellLine,"/frequent_clones.Rds"))
  x2 <- wrap_neighbor_fitness(x,xo,use_gdata = F)
  xo$id <- "frequent"
  x2$id <- "nn"
  xo <- rbind(xo,x2)
  saveRDS(xo,paste0("salehi_data/03_inference_output/",cellLine,"/nn_clones.Rds"))
  
  xmat <- do.call(rbind,lapply(rownames(xo), function(i){
    as.numeric(unlist(strsplit(i,split="[.]")))
  }))
  y <- xo$f_est
  fit <- Krig(xmat,y,m=1)
  saveRDS(fit,file = paste0("salehi_data/03_inference_output/",cellLine,"/krig.Rds"))
}

ff <- list.files("salehi_data/02_optim_data/")
cellLines <- sapply(ff,function(fi) head(unlist(strsplit(fi,split=".Rds")),1))
lapply(cellLines,landscape_inference)

