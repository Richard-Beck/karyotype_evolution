rbf <- function(r,eta=1e-9){
  if(r<eta) return(0)
  return(r^2*log(r))
}

## I think this does the rank reduction as explained in Wood (2003)
eigenpoly <- function(x,f,k){
  d <- as.matrix(dist(x))
  E <- do.call(rbind,lapply(1:nrow(d), function(i){
    sapply(d[i,],rbf)
  }))
  
  Tt <- cbind(1,x)
  
  
  ##potential issue - are we losing negative eigenvalues?
  ## R orders eigenvalues by size but perhaps it was supposed
  ## to be ordered by magnitude. 
  
  ##single value decomposition:
  Uk <- eigen(E)$vectors[,1:k]
  Dk <- (diag(nrow(E))*eigen(E)$values)[1:k,1:k]
  ##Ek:=Uk%*%Dk%*%Uk is the 'optimal' rank k approximation of E. 
  
  ##IDK what happens here
  A0 <- cbind(Uk%*%Dk,Tt)
  B <- t(Tt)%*%Uk
  z <- matrix(0,nrow(B),ncol(Tt))
  
  A <- rbind(A0,cbind(B,z))
  fe <- c(f,rep(0,nrow(B)))
  
  coef <- qr.solve(A,fe)
  
  ## wood says Uk%*%wk~=w where wk is the reduced coefficients.
  ## Therefore, presumably if we extract the approximate to w
  ## we can use this in the existing prediction routine. 
  
  list(coef=coef,w=Uk%*%coef[1:k],v=tail(coef,length(coef)-k),
       knots=x)
}

polypredict <- function(y,ip){
  r <- apply(ip$knots,1,function(i){
    di <- sqrt(sum((i-y)^2))
    rbf(di)
  })
  sum(c(r*ip$w,c(1,y)*ip$v))
}

check_k <- function(k,x,f){
  d <- as.matrix(dist(x))
  E <- do.call(rbind,lapply(1:nrow(d), function(i){
    sapply(d[i,],rbf)
  }))
  
  Tt <- cbind(1,x)
  
  
  ##potential issue - are we losing negative eigenvalues?
  ## R orders eigenvalues by size but perhaps it was supposed
  ## to be ordered by magnitude. 
  
  ##single value decomposition:
  Uk <- eigen(E)$vectors[,1:k]
  Dk <- (diag(nrow(E))*eigen(E)$values)[1:k,1:k]
  ##Ek:=Uk%*%Dk%*%Uk is the 'optimal' rank k approximation of E. 
  
  ##IDK what happens here
  A0 <- cbind(Uk%*%Dk,Tt)
  Ak <- A0%*%solve(t(A0)%*%A0)%*%t(A0)
  B <- t(Tt)%*%Uk
  z <- matrix(0,nrow(B),ncol(Tt))
  A <- rbind(A0,cbind(B,z))
  fe <- c(f,rep(0,nrow(B)))
  coef <- qr.solve(A,fe)
  fp <- A0%*%coef
  
  MSEk <- sum((f-fp)^2)
  tAk <- sum(1-diag(Ak))
  GCVk <- MSEk/(tAk/length(f))^2
  GCVk
}


find_best_k <- function(x,y){
  k_range <- 5:length(y)
  #k_vals <- pbsapply(k_range,tst_k, x=x,f=y,nsplits=nsplits)
  k_vals <- sapply(k_range,function(ki){
    tryCatch(check_k(ki,x=x,f=y),error=function(e) return(Inf))
  })
  k_range[which.min(k_vals)]
}

df <- readRDS("proc_data/optsimple.Rds")
df <- df[df$err<1e3,]
x <- do.call(rbind,lapply(rownames(df), function(r){
  as.numeric(unlist(strsplit(r,split="[.]")))
}))
y <- df$f_est
k <- find_best_k(x,y)
ip <- eigenpoly(x,y,k)
saveRDS(ip,"proc_data/polyscape.Rds")

