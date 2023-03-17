library(transport)
library(ggplot2)
library(pbapply)



setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")

source("rscripts/comparison_functions.R")

f1 <- "ABM/output/v00_rep_005/"
f1 <- paste0(f1,list.files(f1))
f2 <- "ABM/validation/sim_forward/v00_rep_005/"
f2 <- paste0(f2,list.files(f2))
f2 <- paste0(f2,"/",list.files(f2))

get_wdg <- function(g){
  s1 <- lapply(f1,load_sim_dat,g=g)
  s2 <- lapply(f2,load_sim_dat,g=g)
  
  
  d <- do.call(rbind,pblapply(s1, function(i){
    sapply(c(s1,s2), function(j){
      get_dwass(i,j)
    })
  }))
  
  d1 <- d[1:length(s1),1:length(s1)]
  d1 <- mean(d1[upper.tri(d1)])
  d2 <- mean(d[,(1+length(s1)):ncol(d)])
  data.frame(d1=d1,d2=d2,g=g)
}



wdg <- do.call(rbind,lapply(seq(200,2000,200),get_wdg))

x <- reshape2::melt(wdg,id.vars="g")

p <- ggplot(x,aes(x=g,y=value,shape=variable))+
  geom_point()
p
