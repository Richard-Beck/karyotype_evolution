root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
source("rscripts/functions.R")
cpp.out.dir <- "ABM/output/"
opt_batch_file <- paste0(root.dir,"rscripts/simpleopt2.R")

ff <- list.files(cpp.out.dir)
batchname <- "v00"
f1 <- paste0(cpp.out.dir,ff[grepl(batchname,ff)])

f2 <- unlist(lapply(f1,function(fi) {
  paste0(fi,"/",list.files(fi))
  }
))

sapply(f2, function(fi){
  print(fi)
  setwd(root.dir)
  setwd(fi)
  tryCatch(source(opt_batch_file),error=function(e) return(0))
})





