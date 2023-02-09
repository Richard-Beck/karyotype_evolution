library(caret)
library(chebpol)

root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
source("rscripts/functions.R")
source("rscripts/regression_functions.R")
cpp.out.dir <- "ABM/output/"
#proc_batch_file <- paste0(root.dir,"rscripts/predict_beyond.R")

ff <- list.files(cpp.out.dir)
batchname <- "secondTest"
f1 <- paste0(cpp.out.dir,ff[grepl(batchname,ff)])

f2 <- unlist(lapply(f1,function(fi) {
  paste0(fi,"/",list.files(fi),"/")
}
))

x <- do.call(rbind,lapply(f2, function(dir){
  fi <- list.files(paste0(dir,"proc_data"))
  if("landscape_metrics.Rds"%in%fi){
    xi <- readRDS(paste0(dir,"proc_data/landscape_metrics.Rds"))
    ids <- unlist(strsplit(dir,split="/"))
    landscape_id <- unlist(strsplit(ids[3],split="_"))[3]
    xi$landscape_id <- landscape_id
    xi$rep_id <- ids[4]
    return(xi)
    
  }
  return(NULL)
}))

p <- ggplot(x,aes(x=ids,y=value,color=method,group=method))+
  facet_grid(rows=vars(variable),cols=vars(landscape_id),scales="free")+
  geom_smooth()
p



