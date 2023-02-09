### GAM example using mgcv

library(mgcv)
library(ggplot2)

gaussf <- function(xi) exp(-sqrt(sum((xi-c(2,2))^2))) + exp(-sqrt(sum((xi-c(7,7))^2))/3^2) 

x <- data.frame(matrix(runif(800,min=0,max=10),ncol=2))
colnames(x) <- c("x0","x1")
y <- apply(x,1, gaussf)

dat <- cbind(x,y=y)

k <- 80
m <- 2
b1 <- mgcv::gam(y ~ s(x0, x1,bs='tp'), data = dat)
summary(b1)
plot(b1)

tst <- expand.grid(seq(0,10,0.5),seq(0,10,0.5))
colnames(tst) <- colnames(x)

tst$z <- apply(tst,1,gaussf)
tst$y <- predict(b1,tst)


tst <- reshape2::melt(tst,id.vars=c("x0","x1"))

p <- ggplot(tst,aes(x=x0,y=x1,fill=value))+
  facet_grid(rows=vars(variable))+
  scale_fill_viridis_c()+
  geom_raster()
p
