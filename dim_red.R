points <- expand.grid(1:10,1:10,1:10)
centre <- c(5,5,5)

d <- apply(points,1, function(p) sqrt(sum((centre-p)^2)))
points$dist <- d
points <- points[d<3,]

library(Rtsne)
library(ggplot2)

x <-Rtsne(points)
x <- data.frame(x$Y)
x$dist <- points$dist
x$x <- points$Var1


p <- ggplot(x,aes(x=X1,y=X2,color=dist))+
  geom_point()+
  scale_color_viridis_c()
p