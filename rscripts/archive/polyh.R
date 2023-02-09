## code and functions for polyharmonic interpolation, implementation based
## off https://en.wikipedia.org/wiki/Polyharmonic_spline

## radial basis functions:
rbf <- function(r,k){
  #return(exp(-(0.05*r)^2))
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
  }))
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

#library(chebpol)
## check functions work and give same answer as chebpol:
#centres <- cbind(c(6,3,1,6),c(5,1,2,6))
#f <- c(0.7,0.2,0.9,0.1)
#k <- 2
#ph <- setup_polyh(centres,f,k)
#ph2 <- ipol(val = f,knots = t(centres),k=2,method = "polyh")
#ph
#poly_interp(c(5,5),centres,ph$w,ph$v)
#ph2(c(5,5))

