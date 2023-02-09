setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")
source("rscripts/polyh.R")
landscape <- read.table("ABM/landscapes/randomTest_rep_001.txt",sep=",")
config <- readLines("ABM/config/randomTest_rep_001.txt")
scale <- config[grepl("scale",config)]
scale <- as.numeric(unlist(strsplit(scale,split="scale,")))[2]


getf <- function(vec){
  d <- apply(landscape,1,function(ci) {
    sqrt(sum((vec-ci)^2))
  })
  sum(sin(d)*scale)
}

nchrom <- ncol(landscape)
library(pbapply)
peaks <- do.call(rbind,pblapply(1:300, function(i){
  init <- sample(1:8,nchrom,replace = T)
  opt <- optim(init,fn=getf,lower=1,method="L-BFGS-B")
  peak <- round(opt$par)
  peak
}))

f <- -apply(peaks,1,getf)
peaks <- apply(peaks,1,paste,collapse=".")

df <- data.frame(peak=peaks,f=f,n=1)
df <- aggregate(list(n=df$n),by=list(peak=df$peak,f=df$f),sum)
