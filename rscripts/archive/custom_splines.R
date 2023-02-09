## One way of constructing a GRF is by assuming that the field is the sum of a 
## large number of plane, cylindrical or spherical waves with uniformly
## distributed random phase. Where applicable, the central limit theorem 
## dictates that at any point, the sum of these individual plane-wave 
##contributions will exhibit a Gaussian distribution. 
library(ggplot2)



create_var <- function(np=20,minv=0,maxv=10) runif(np,minv,maxv)
create_vars <- function(nvar=1,np=20,minv=0,maxv=10){
  centres <- data.frame(do.call(cbind,lapply(1:nvar,function(dummy){
    create_var(np=np,minv=minv,maxv=maxv)
  })))
  colnames(centres) <- paste0("x",1:ncol(centres))
  return(centres)
}

getf <- function(var,centres,noise=0){
  d <- apply(centres,1,function(ci) {
    sqrt(sum((var-ci)^2))
  })
  f <- sin(d*freq)+rnorm(1,mean=0,sd=noise)
  sum(f)
}


rbf <- function(d){
  d <- max(d,1e-20)
  return(d^0.5)
  if(d==0) return(d^2)
  d^2*log(d)
}

mrb <- function(x1,x2,C){
  # -0.9080199
  s0 <- sapply(1:nrow(x1),function(i){
    sum(sapply(1:nrow(x2), function(j){
      dij <- sum((x1[i,]-x2[j,])^2)
      rbf(dij)
    })*C)
  })
  return(s0)
  
}

## should return the same as the tps function
rich.predict <- function(opt,y){
  object <- opt
  xx <- y
  xc <- object$transform$x.center
  xs <- object$transform$x.scale
  
  #tst <- scale(xx,xc,xs)
  xx <- do.call(cbind,lapply(1:ncol(xx), function(i){
    (xx[,i]-xc[i])/xs[i]
  }))
  #xx <- scale(xx, xc, xs)
  temp.c <- object$c
  temp.d <- object$d
  Tmatrix <- cbind(1,xx)
  temp <- Tmatrix %*% temp.d[object$ind.drift]
  d <- ncol(xx)
  m <- (d+1)/2
  
  c1 <- mrb(x1 = xx, x2 = object$knots, C = temp.c)*radbas.constant(m,d)
  c2 <- do.call(object$cov.function.name, c(object$args, 
                                            list(x1 = xx, x2 = object$knots, C = temp.c)))
  #46.404591  -6.349416  -2.816273 -27.089015  -6.396628  55.1909740
  temp <- temp + c1
  return(temp)
}

set.seed(42)



args <- list(p=1)
x1 <- matrix(1:5,nrow=5)
x2 <- x1
C <- x1
d <- 1
n1 <- 5
n2 <- 5
n3 <- 1
par <- c(0.5,0)

y0 <- .Fortran("multrb", PACKAGE = "fields", nd = as.integer(d), 
                 x1 = as.double(x1), n1 = as.integer(n1), x2 = as.double(x2), 
                 n2 = as.integer(n2), par = as.double(par), c = as.double(C), 
                 n3 = as.integer(n3), h = as.double(rep(0, n1 * 
                                                          n3)), work = as.double(rep(0, n2)))$h

y1 <- do.call("Rad.cov",c(args, 
                    list(x1 = x1, x2 = x2, C = C)))

m <- (d+1)/2

y2 <- mrb(x1,x2,C)*radbas.constant(m,d)

nvar <- 3
centres <- create_vars(nvar=nvar,np=5)
x <- create_vars(nvar=nvar,np=100,minv=2,maxv=8)
fx <- apply(x,1,getf,centres=centres,noise=0.1)

library(fields)
opt <- Tps(x,fx)
fpred <- predict(opt,x)
fpred2 <- rich.predict(opt,x)
plot(fpred,fx)
plot(fpred,fpred2)

y <- data.frame(x1=seq(0,10,0.1))
f_tru <-  apply(y,1,getf,centres=centres,noise=0)
f_est <-  predict(opt,y)
f_est2 <- rich.predict(opt,y)

y$f_tru <- f_tru
y$f_est <- f_est

x$f <- fx

y <- reshape2::melt(y,measure.vars=c("f_tru","f_est"))


p <- ggplot(x,aes(x=x1))+
  geom_point(aes(y=f))+
  geom_line(data=y,aes(y=value,color=variable,group=variable))+
  scale_x_continuous("chromosome 1")+
  scale_y_continuous("fitness")
p
