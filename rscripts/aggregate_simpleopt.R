root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
source("rscripts/functions.R")
cpp.out.dir <- "ABM/output/"

ff <- list.files(cpp.out.dir)
batchname <- "randomTest"
f1 <- paste0(cpp.out.dir,ff[grepl(batchname,ff)])

f2 <- unlist(lapply(f1,function(fi) {
  paste0(fi,"/",list.files(fi))
  }
))

x <- lapply(1:length(f2), function(i){
  fi <- f2[i]
  print(fi)
  setwd(root.dir)
  setwd(fi)
  if(!"optsimple.Rds"%in%list.files("proc_data")) return(NULL)
  df <- readRDS("proc_data/optsimple.Rds")
  df$id <- i
  return(df)
})

df <- do.call(rbind,x)
library(ggplot2)

p <- ggplot(df,aes(x=f_est,y=f_tru,color=err))+
  scale_color_viridis_c(trans="log")+
  geom_point()+
  geom_abline(color="red")+
  scale_x_continuous("estimated fitness")+
  scale_y_continuous("true fitness")
p

p <- ggplot(df[df$err<10^6,],aes(x=f_est,y=f_tru))+
  geom_point()+
  geom_abline(color="red",size=2)+
  scale_x_continuous("estimated fitness")+
  scale_y_continuous("true fitness")
p



