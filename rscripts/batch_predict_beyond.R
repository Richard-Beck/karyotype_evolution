library(caret)
library(chebpol)

root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
source("rscripts/functions.R")
source("rscripts/regression_functions.R")
cpp.out.dir <- "ABM/output/"
proc_batch_file <- paste0(root.dir,"rscripts/predict_beyond.R")

ff <- list.files(cpp.out.dir)
batchname <- "secondTest"
f1 <- paste0(cpp.out.dir,ff[grepl(batchname,ff)])

f2 <- unlist(lapply(f1,function(fi) {
  paste0(fi,"/",list.files(fi),"/")
}
))

sapply(f2, function(dir){
  setwd(root.dir)
  print(dir)
  setwd(dir)
  tryCatch(source(proc_batch_file),error=function(e) print("error!"))
})

