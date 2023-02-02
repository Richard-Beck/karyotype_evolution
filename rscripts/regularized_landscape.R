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
setup_polyh <- function(x,B,lambda){
  M <- as.matrix(dist(x))
  M <- do.call(rbind,lapply(1:nrow(M), function(i){
    sapply(1:ncol(M), function(j) rbf(M[i,j],k))
  }))
  lambda_mat <- diag(nrow(M))*lambda
  solve(t(M)%*%M+lambda)%*%t(M)%*%B
}


## interpolate a new point
poly_interp <- function(pt,x,theta){
  M <- apply(knots,1, function(ki) rbf(sqrt(sum((ki-pt)^2)),k))
  M%*%theta
}

## write out interpolator to use in cpp
write_poly <- function(path,knots,w,v){
  out <- cbind(knots,w)
  out <- rbind(out,v)
  options(scipen=999)
  write.table(x = out,file = path,sep=",",col.names = F,row.names = F)
}

library(chebpol)

distf <- function(a,b) sqrt(sum((a-b)^2))

## check functions work and give same answer as chebpol:
npoints <- 4
centres <- cbind(runif(npoints),runif(npoints))
f <- apply(centres,1,function(i) exp(-distf(i,c(0.2,0.2))))
ph <- setup_polyh(centres,f,lambda=0)
apply(centres,1,poly_interp)




