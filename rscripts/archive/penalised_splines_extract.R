library(mgcv)

distf <- function(a,b) sqrt(sum((a-b)^2))
## check functions work and give same answer as chebpol:
npoints <- 50
centres <- cbind(runif(npoints),runif(npoints))
f <- apply(centres,1,function(i) exp(-distf(i,c(0.2,0.2))))
dat <- data.frame(cbind(centres,f))


b <- gam(f ~ s(V1,V2,bs='tp'), data = dat)

newd <- data.frame(V1=0.33,V2=0.53)
h <- predict(b, newd)
#Xf <- pr3(b$smooth[[1]], as.list(newd))
#Xf <- t(qr.qty(qrc, t(X))[(j + 1):k, , drop = FALSE])
Xf <- PredictMat(b$smooth[[1]],newd)
X <- c(1,Xf)
X%*%b$coefficients
