setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")
source("rscripts/polyh.R")
landscape <- read.table("ABM/landscapes/randomTest_rep_001.txt",sep=",")
config <- readLines("ABM/config/randomTest_rep_001.txt")
opt <- readRDS("ABM/output/randomTest_rep_001/00000/proc_data/optsimple.Rds")
opt <- opt[opt$err<1e3,]
opt <- opt[order(opt$f_est,decreasing=T),]
centres <- do.call(rbind,lapply(rownames(opt), function(oi){
  pts <- as.numeric(unlist(strsplit(oi,split="[.]")))
  ##sometimes odd stuff happens (we could get a column of all ones)
  ## is adding a little bit of noise a reasonable solution?
  pts <- pts + runif(length(pts),min=-0.02,max=0.02)
}))
f <- opt$f_est

ph <- setup_polyh(knots=centres,values=f,k=2)


peak <- centres[2,]

scale <- config[grepl("scale",config)]
scale <- as.numeric(unlist(strsplit(scale,split="scale,")))[2]

#examine <- "peak"


for(i in 1:9){
  grid <- expand.grid(x=seq(0,8,1),y=seq(0,8,1))
  chr <- c(i,i+1)
  f <- sapply(1:nrow(grid),function(i){
    vec <- peak
    if(examine == "diploid"){
      vec <- rep(2,length(vec))
    }
    vec[chr] <- as.numeric(grid[i,])
    d <- apply(landscape,1,function(ci) {
      sqrt(sum((vec-ci)^2))
    })
    f_pred <- poly_interp(vec,centres,ph$w,ph$v)
    return(c(f_pred=f_pred,f_tru=sum(sin(d)*scale)))
  })
  f <- data.frame(t(f))
  grid <- cbind(grid,f)
  
  grid <- reshape2::melt(grid,id.vars=c("x","y"))
  
  library(ggplot2)
  
  chrx <- peak[-chr]
  centrex <- centres[,-chr]
  hits <- apply(centrex,1, function(xx) mean(xx==chrx)==1)
  hits <- data.frame(centres)[hits,chr]
  colnames(hits) <- c("x","y")
  
  p <- ggplot(grid,aes(x=x,y=y))+
    facet_grid(cols=vars(variable))+
    geom_raster(aes(fill=value))+
    geom_point(data=hits)+
    scale_fill_viridis_c()+
    scale_x_continuous(paste("chromosome",chr[1]),breaks=0:8)+
    scale_y_continuous(paste("chromosome",chr[2]),breaks=0:8)
  plot(p)
}






