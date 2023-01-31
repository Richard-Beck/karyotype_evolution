pk1 <- c(3,3)
pk2 <- c(7,7)
sigma <- 2
gaussf <- function(xy){

  max(exp(-sum((xy-pk1)^2/sigma^2)),
      exp(-sum((xy-pk2)^2/sigma^2)))
}



trn <- do.call(rbind,lapply(1:50, function(nothing){
  xy <- runif(2,min=0,max=10)
  f <- gaussf(xy)
  data.frame(x=xy[1],y=xy[2],f)
}))

library(chebpol)
tst <- expand.grid(1:10,1:10)
colnames(tst) <- c("x","y")
s <- seq(0,1,length.out=10)

grid <- list(s,s)
colnames(tst) <- c("x","y")
knots <- t(as.matrix(trn[,1:2]))
ph2 <- ipol(trn[,3], grid=grid,knots=knots,  method='polyh')
m <- t(as.matrix(tst[,1:2]))
m <- apply(m,2,as.numeric)

tst$f <- ph2(m)

wrap_ph <- function(pars){
  d1 <- 0
  d2 <- 0
  if(max(pars)>10) d1 <- max(pars)-10
  if(min(pars)<0) d2 <- 0-min(pars)
  return(-ph2(pars)+d1+d2)
} 



p <- ggplot(trn,aes(x=x,y=y,color=f))+
  geom_point()+
  scale_color_viridis_c()
p
p <- ggplot(tst,aes(x=x,y=y,fill=f))+
  geom_raster()+
  scale_fill_viridis_c()
p

maxima <- data.frame(do.call(rbind,lapply(1:100, function(i){
  p0 <- runif(2,0,10)
  optim(p0,wrap_ph)$par
})))
colnames(maxima)=c("x","y")

p <- ggplot(maxima,aes(x=x,y=y))+
  geom_jitter(height=0.3,width=0.3)+
  scale_x_continuous(limits=c(0,10))+
  scale_y_continuous(limits=c(0,10))
p

