## One way of constructing a GRF is by assuming that the field is the sum of a 
## large number of plane, cylindrical or spherical waves with uniformly
## distributed random phase. Where applicable, the central limit theorem 
## dictates that at any point, the sum of these individual plane-wave 
##contributions will exhibit a Gaussian distribution. 
library(ggplot2)
type <- "spherical"
if(type=="spherical"){
  
  centres <- data.frame(x=runif(10,0,10),y=runif(10,0,10))
  freq <- pi/2
  
  grid <- expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))
  
  f <- do.call(rbind,lapply(1:nrow(grid),function(i){
    d <- apply(centres,1,function(ci) sqrt(sum((grid[i,]-ci)^2)))
    f <- sin(d*freq)
    data.frame(f1=f[1],f2=f[2],sum.10=sum(f))
  }))
  
  grid <- cbind(grid,f)
  grid <- reshape2::melt(grid,id.vars=c("x","y"))
  
  p <- ggplot(grid,aes(x=x,y=y,fill=value))+
    facet_grid(cols=vars(variable))+
    scale_fill_viridis_c("fitness")+
    scale_x_continuous("chromosome 1")+
    scale_y_continuous("chromosome 2")+
    geom_raster()
  p
  
  ## an interesting is what we expect the fitness distribution to be for 
  ## N waves
  N <- seq(10,300,10)
  maxf <- sapply(N, function(n){
    max(abs(sapply(1:10000, function(i) sum(sin(runif(n,min=0,max=2*pi))))))
    
  })
  
  plot(N,maxf)
  lines(N,pi*sqrt(N))
  
  
}





