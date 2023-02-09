
centres <- data.frame(x=runif(5,0,10),y=runif(5,0,10))
freq <- pi/8

ff <- function(xy){
  d <- apply(centres,1,function(ci) sqrt(sum((xy-ci)^2)))
  f <- sin(d*freq)
  sum(f)
}


trn <- do.call(rbind,lapply(1:10, function(nothing){
  xy <- runif(2,min=2,max=8)
  f <- ff(xy)
  data.frame(x=xy[1],y=xy[2],f)
}))

library(chebpol)
tst <- expand.grid(1:10,1:10)
colnames(tst) <- c("x","y")
colnames(tst) <- c("x","y")
knots <- t(as.matrix(trn[,1:2]))
ph2 <- ipol(trn[,3], knots=knots,  method='polyh')
m <- t(as.matrix(tst[,1:2]))
m <- apply(m,2,as.numeric)
tst$f <- ph2(m)
tst$tru <- apply(tst[,1:2],1,ff)
tst <- reshape2::melt(tst,id.vars=c("x","y"))

p <- ggplot(tst,aes(x=x,y=y))+
  facet_grid(cols=vars(variable))+
  geom_raster(aes(fill=value))+
  geom_point(data=trn)+
  scale_fill_viridis_c()
p

