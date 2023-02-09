setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")
source("rscripts/polyh.R")
landscape <- read.table("ABM/landscapes/randomTest_rep_001.txt",sep=",")
config <- readLines("ABM/config/randomTest_rep_001.txt")

scale <- config[grepl("scale",config)]
scale <- as.numeric(unlist(strsplit(scale,split="scale,")))[2]

opt <- readRDS("ABM/output/randomTest_rep_001/00000/proc_data/optsimple.Rds")
opt <- opt[opt$err<1e3,]
opt <- opt[order(opt$f_est,decreasing=T),]
centres <- rownames(opt)

getf <- function(vec){
  d <- apply(landscape,1,function(ci) {
    sqrt(sum((vec-ci)^2))
  })
  sum(sin(d)*scale)
}

gen_all_neighbours <- function(ki){
  unlist(lapply(1:length(ki), function(i){
    ki1 <- ki
    ki2 <- ki
    ki1[i] <- k[i]+1
    ki2[i] <- k[i]-1
    c(paste(ki1,collapse="."),paste(ki2,collapse="."))
  }))
}

k <- do.call(rbind,lapply(centres, function(ci){
  as.numeric(unlist(strsplit(ci,split="[.]")))
}))

df <- do.call(rbind,lapply(1:nrow(k),function(i){
  x <- gen_all_neighbours(k[i,])
  x <- x[!x%in%rownames(opt)]
  xn <- lapply(x, function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  })
  f_tru <- sapply(xn,getf)
  df <- do.call(rbind,lapply(1:length(x), function(dummy) opt[i,]))
  rownames(df) <- x
  df$f_est <- df$f_est*0.99
  df$f_tru <- f_tru
  df
}))
df <- df[!duplicated(rownames(df)),]

plot(opt$f_est,opt$f_tru)
plot(df$f_est,df$f_tru)

df <- rbind(opt,df)
saveRDS(df,"ABM/output/randomTest_rep_001/00000/proc_data/optsimpletst.Rds")
