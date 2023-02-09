## code and functions for polyharmonic interpolation, implementation based
## off https://en.wikipedia.org/wiki/Polyharmonic_spline

## radial basis functions:
rbf <- function(r,k){
  val <- NaN
  if(r>=1){
    if(k%%2==0){
      val <- r^k*log(r)
    }else{
      val <- r^k
    }
  }else{
    if(k%%2==0){
      val <- r^(k-1)*log(r^r)
    }else{
      val <- r^k
    }
  }
  return(val)
}

##setup the interpolator, compute w and v
setup_polyh <- function(knots,values,k){
  d <- as.matrix(dist(knots))
  A <- do.call(rbind,lapply(1:nrow(d), function(i){
    sapply(1:ncol(d), function(j) rbf(d[i,j],k))
  })) ## this is the design matrix
  
  B <- cbind(1,knots)
  
  z0 <- matrix(0,ncol(B),ncol(B))
  
  m <- rbind(cbind(A,B),cbind(t(B),z0))
  f <- c(values,rep(0,ncol(B)))
  wv <- solve(m,f)
  
  w <- wv[1:nrow(knots)]
  v <- wv[(nrow(knots)+1):length(wv)]
  return(list(w=w,v=v))
}


## interpolate a new point
poly_interp <- function(pt,knots,w,v,k=2){
  p2 <- c(1,pt)
  r <- apply(knots,1, function(ki) rbf(sqrt(sum((ki-pt)^2)),k))
  sum(w*r)+sum(v*p2)
}

## write out interpolator to use in cpp
write_poly <- function(path,knots,w,v){
  out <- cbind(knots,w)
  out <- rbind(out,v)
  options(scipen=999)
  write.table(x = out,file = path,sep=",",col.names = F,row.names = F)
}

test_subset <- function(n,centres,f,k=2){
  keep <- sample(1:length(f),n)
  tst <- (1:length(f))[-keep]
  Xtrain <- centres[keep,]
  Ytrain <- f[keep]
  Xtest <- centres[tst,]
  Ytest <- f[tst]
  ph <- setup_polyh(Xtrain,Ytrain,k)
  Ypred <- apply(Xtest,1,poly_interp,knots=Xtrain,w=ph$w,v=ph$v,k)
  errs <- sqrt(mean((Ypred-Ytest)^2))
  cbind(keep,errs)
}

remove_1_point <- function(centres,f,k=2){
  
  errs <- sapply(1:nrow(centres), function(i){
    cx <- centres[-i,]
    fx <- f[-i]
    ph <- setup_polyh(cx,fx,k)
    (poly_interp(centres[i,],cx,ph$w,ph$v,k)-f[i])^2
  })
  
  cx <- centres[-which.min(errs),]
  fx <- f[-which.min(errs)]
  list(centres=cx,f=fx)
  
}

library(ggplot2)

distf <- function(a,b) sqrt(sum((a-b)^2))
dxy <- 0.1
lscape <- expand.grid(x=seq(0,1,dxy),y=seq(0,1,dxy))
lscape$f <- apply(lscape,1,function(i) exp(-distf(i,c(0.2,0.2))))

## check functions work and give same answer as chebpol:
npoints <- 30
n_knots <- 20
centres <- cbind(runif(npoints),runif(npoints))
f <- apply(centres,1,function(i) exp(-distf(i,c(0.2,0.2))))
f <- f+rnorm(length(f),mean=0,sd=0.04)
k <- 2

errs <- rep(NaN,length(f)-2)

ph <- setup_polyh(centres,f,k=2)
e0 <- sqrt(mean((sapply(1:nrow(centres), function(i) poly_interp(centres[i,],centres,ph$w,ph$v,k=2))-f)^2))

errs[1] <- e0

ds <- list(centres=centres,f=f)

for(i in 2:length(errs)){
  ds <- remove_1_point(ds$centres,ds$f,k=2)
  ph2 <- setup_polyh(ds$centres,ds$f,k=2)
  errs[i] <- sqrt(mean((sapply(1:nrow(centres), function(i) poly_interp(centres[i,],ds$centres,ph2$w,ph2$v,k=2))-f)^2))
  
}

plot(errs)
ds <- list(centres=centres,f=f)
for(i in 1:20){
  ds <- remove_1_point(ds$centres,ds$f,k=2)
}

ph2 <- setup_polyh(ds$centres,ds$f,k=2)
lscape$f_est <- apply(lscape[,1:2],1,function(i) poly_interp(i,ds$centres,ph2$w,ph2$v,k=2))
lscape$f_est_full <- apply(lscape[,1:2],1,function(i) poly_interp(i,centres,ph$w,ph$v,k=2))
lscape <- reshape2::melt(lscape,id.vars=c("x","y"))

df1 <- data.frame(centres)
df2 <- data.frame(ds$centres)
colnames(df1) <- c("x","y")
df1$variable <- "f_est_full"
colnames(df2) <- c("x","y")
df2$variable <- "f_est"
p <- ggplot(lscape,aes(x=x,y=y))+
  facet_grid(cols=vars(variable))+
  geom_raster(aes(fill=value))+
  geom_point(data=df1)+
  geom_point(data=df2)+
  scale_fill_viridis_c()
p

if(FALSE){
  
  trn <- data.frame(do.call(rbind,lapply(1:100, function(i) test_subset(n_knots,centres,f,k=2))))
  
  errs <- aggregate(list(errs=trn$errs),by=list(keep=trn$keep),mean)
  keep <- (errs$keep[order(errs$err)])[1:n_knots]
  
  df <- data.frame(centres)
  colnames(df) <- c("x","y")
  
  df$kept <- "false"
  df$kept[keep] <- "true"
  
  x <- centres[keep,]
  y <- f[keep]
  
  ph <- setup_polyh(x,y,k=2)
  ph_full <- setup_polyh(centres,f,k=2)
  
  lscape$f_est <- apply(lscape[,1:2],1,function(i) poly_interp(i,x,ph$w,ph$v,k=2))
  lscape$f_est_full <- apply(lscape[,1:2],1,function(i) poly_interp(i,centres,ph_full$w,ph_full$v,k=2))
  lscape <- reshape2::melt(lscape,id.vars=c("x","y"))
  
  p <- ggplot(lscape,aes(x=x,y=y))+
    facet_grid(cols=vars(variable))+
    geom_raster(aes(fill=value))+
    geom_point(data=df,aes(color=kept))+
    scale_fill_viridis_c()
  p
  
}







