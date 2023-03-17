##prepare data to fit sim
proc_sim <- function(dir,times){
  dt <- read.csv(paste0(dir,"/log.txt"),sep=",")$dt
  
  filenames <- list.files(dir)
  filenames <- filenames[!filenames%in%c("proc_data","log.txt")]
  tx <- as.numeric(substr(filenames,1,5))
  filenames <- paste0(dir,"/",sapply(times, function(ti) filenames[which.min(abs(ti-tx))]))
  tx <- sapply(times, function(ti) tx[which.min(abs(ti-tx))])
  
  x <- lapply(filenames,function(fi){
    xi <- read.csv(fi,header=F)
    nchrom <- ncol(xi)-2
    k <- apply(xi[,1:nchrom],1,paste,collapse=".")
    n <- xi[,nchrom+1]
    fitness <- xi[,nchrom+2]
    data.frame(karyotype=k,n=n,fitness=fitness)
  })
  
  
  
  pop.fitness <- sapply(x,function(xi) sum(xi$n*xi$fitness)/sum(xi$n))
  
  unique.karyotypes <- unique(do.call(rbind,x)$karyotype)
  
  x <- do.call(rbind,lapply(unique.karyotypes, function(k){
    sapply(x, function(xi) sum(xi$n[xi$karyotype==k]))
  }))
  
  colnames(x) <- tx
  rownames(x) <- unique.karyotypes  
  list(x=x,pop.fitness=pop.fitness,dt=dt)
}

optimx <- function(f,u,fit_dat){
  
  ux <- u/colSums(fit_dat$x)
  up <- 1-ux
  #ux <- pmax(0.001,pmin(ux,0.999))
  #up <- pmax(0.001,pmin(up,0.999))
  ux0 <- ux
  up0 <- up
  fp <- (fit_dat$pop.fitness-ux*f)/up#pmin(3,pmax(0,(fit_dat$pop.fitness-ux*f)/up))
  
  tx <- as.numeric(colnames(fit_dat$x))
  
  ## we would like to capture any fitness implications of the first observable timepoint
  i <- min(which(u>0))
  if(FALSE & i>1 & length(which(u>0))==1){
    ## assume that a clone emerged in the middle of the observation period"
    t <- fit_dat$dt*diff(tx)[i-1]/2
    ux[i] <- ux0[i]/1000*exp(t*f)
    up[i] <- max(0.0001,up0[i-1])*exp(t*mean(c(fp[i],fp[i-1])))
    sx <- ux[i]+up[i]
    ux[i]<- ux[i]/sx
    up[i]<- up[i]/sx
  }
  
  for(i in (1+min(which(u>0))):length(ux)){
    t <- fit_dat$dt*diff(tx)[i-1]
    ux[i] <- ux0[i-1]*exp(t*f)
    up[i] <- max(0.0001,up0[i-1])*exp(t*mean(c(fp[i],fp[i-1])))
    sx <- ux[i]+up[i]
    ux[i]<- ux[i]/sx
    up[i]<- up[i]/sx
  }
  #print(ux)
  negll <- -log(dbinom(u,colSums(fit_dat$x),prob=ux))
  negll[!is.finite(negll)] <- 10^9
  sum(negll)
}

opt_tp1 <- function(ftp1,n1,x1,f1,ff,t){
  ux0 <- x1/sum(c(x1,n1))
  up0 <- 1-ux0
  
  fp <- (ff-ux0*ftp1)/up0
  
  ux <- ux0*exp(t*ftp1)
  up <- pmax(0.0001,up0)*exp(t*fp)
  sx <- ux+up
  ux <- ux/sx
  
  negll <- -log(dbinom(0,sum(c(n1,x1)),prob=ux))
  f_all <- c(ftp1,f1)
  negll[!is.finite(negll)] <- 10^9
  negd <- -log(dnorm(ftp1,mean(f_all),max(0.001,sd(f_all))))
  sum(negll)+sum(negd)
}

optimsimple_old <- function(x){
  ntp <- apply(x$x,1,function(i) sum(i>0))
  n <- rowSums(x$x)
  to_fit <- which(n>5)
  x$x <- x$x[to_fit,,drop=FALSE]
  ##optimx function doesn't work when 
  ##hmm..
  ##dompoints <- apply(x$x,2,function(xi) sum(xi>0)==1)
  if(length(to_fit)==1){
    x_opt <- data.frame(f_est=mean(x$pop.fitness),err=0,n=sum(x$x),ntp=ncol(x$x))
    rownames(x_opt) <- rownames(x$x)
    return(x_opt)
  }
  ix <- 1:nrow(x$x)
  
  peakest <- apply(x$x,1,which.max)
  init_guess <- x$pop.fitness[peakest]
  opt <- do.call(rbind,lapply(ix, function(i){
    u <- x$x[i,]
    fi <- u/colSums(x$x)
    fmax <- max(fi)
    opt <- NULL
    if(fmax>0.95){
      ## if there is only one clone present the optimx function breaks,
      ## so it is better to just use the pop growth rate to determine
      ## fitness of that clone
      opt <- list(minimum=x$pop.fitness[which.max(fi)],objective=0)
    }else{
      igi <- init_guess[i]
      val <- Inf
      ntries <- 1
      
      while(val>1e4 & ntries < 10){
        ntries <- ntries + 1
        interval <- c(igi-1/ntries,igi+1/ntries)
        opt <- optimise(optimx,interval=interval,u=u,fit_dat=x)
        val <- opt$objective
      }
      if(opt$objective>1e6) opt$minimum <- igi
    }
    
    data.frame(f_est=opt$minimum,err=opt$objective,
               n=sum(u),ntp=sum(u>0))
    
  }))
  rownames(opt) <- rownames(x$x)
  opt
  #  opt[opt$err<1e4,]
  
}

optimsimple <- function(x){
  ntp <- apply(x$x,1,function(i) sum(i>0))
  tp1 <- x$x[,1]>0 & ntp==1
  n <- rowSums(x$x)
  to_fit <- which(n>5 & !tp1)
  xtp1 <- x$x[n>5 & tp1,,drop=FALSE]
  x$x <- x$x[to_fit,,drop=FALSE]
  
  
  if(length(to_fit)==1){
    x_opt <- data.frame(f_est=mean(x$pop.fitness),err=0,n=sum(x$x),ntp=ncol(x$x))
    rownames(x_opt) <- rownames(x$x)
    return(x_opt)
  }
  ix <- 1:nrow(x$x)
  
  peakest <- apply(x$x,1,which.max)
  init_guess <- x$pop.fitness[peakest]
  opt <- do.call(rbind,lapply(ix, function(i){
    u <- x$x[i,]
    fi <- u/colSums(x$x)
    fmax <- max(fi)
    opt <- NULL
    if(fmax>0.95){
      ## if there is only one clone present the optimx function breaks,
      ## so it is better to just use the pop growth rate to determine
      ## fitness of that clone
      opt <- list(minimum=x$pop.fitness[which.max(fi)],objective=0)
    }else{
      igi <- init_guess[i]
      val <- Inf
      ntries <- 1
      
      while(val>1e4 & ntries < 10){
        ntries <- ntries + 1
        interval <- c(igi-1/ntries,igi+1/ntries)
        opt <- optimise(optimx,interval=interval,u=u,fit_dat=x)
        val <- opt$objective
      }
      if(opt$objective>1e6) opt$minimum <- igi
    }
    
    data.frame(f_est=opt$minimum,err=opt$objective,
               n=sum(u),ntp=sum(u>0))
    
  }))
  rownames(opt) <- rownames(x$x)
  if(nrow(xtp1)==0) return(opt)
  ff <- mean(x$pop.fitness[1:2])
  tp1 <- x$x[,1]>0 
  n1 <- x$x[,1]
  x1 <- xtp1[,1]
  f1 <- opt$f_est[n1>0]
  ftp1 <- rep(x$pop.fitness[1],nrow(xtp1))
  t <- fit_dat$dt*diff(as.numeric(colnames(xtp1)[1:2]))
  opt1 <- data.frame()
  if(length(x1)==1){
    opt1 <- optimise(opt_tp1,lower=-1,upper=ff,n1=n1,x1=x1,f1=f1,ff=ff,t=t)
    opt1 <- data.frame(f_est=opt1$minimum,err=opt1$objective,n=x1,ntp=1)
  }else{
    opt1 <- optim(ftp1,fn=opt_tp1,n1=n1,x1=x1,f1=f1,ff=ff,t=t)
    opt1 <- data.frame(f_est=opt1$par,err=opt1$value,n=x1,ntp=1)
    
  }
  rownames(opt1) <- rownames(xtp1)
  opt <- rbind(opt,opt1)
  opt
}

## function for inferring nn fitness
optimf_nn <- function(par,u,tmat,known_rates){
  
  est_rates <- par[1:nrow(tmat)]
  est_rates <- pmin(est_rates,0.99*min(known_rates)) ## problems if multiple known rates
  distr_par <- tail(par,2)
  ll1 <- -sum(dnorm(c(est_rates,known_rates),distr_par[1],distr_par[2],log=T))
  
  rmat <- tmat%*%u[colnames(tmat),]
  ll2 <- sum(sapply(1:nrow(rmat), function(i){
    idi <- rownames(tmat)[i]
    ## this now fails when there are multiple columns in tmat:
    upred <- rmat[idi,]*known_rates/(known_rates-est_rates[i])
    u_obs <- rep(0,length(upred))
    if(idi%in%rownames(u)) u_obs <- u[idi,]
    -sum(dpois(u_obs,upred,log=T))
  }))
  
  ll1 + ll2
  
}

infer_nn <- function(x,x_opt,pm0){
  u <- x$x
  nx <- do.call(rbind,lapply(rownames(x_opt), function(i) as.numeric(unlist(strsplit(i,split="[.]")))))
  nn <- gen_all_neighbours(rownames(x_opt))
  tmat <- apply(nx,1,function(xi){
    apply(nn,1,function(xj){
      prod(sapply(1:length(xj), function(k) pij(xi[k],xj[k],pm0)))
    })
  })
  
  colnames(tmat) <- rownames(x_opt)
  rownames(tmat) <- apply(nn,1,paste,collapse=".")
  known_rates <- mean(x_opt$f_est)
  par <- c(rep(0,nrow(tmat)),mean(x_opt$f_est),1)
  opt <- optim(par,optimf_nn,u=u,tmat=tmat,known_rates=known_rates)
  
  dfstat <- do.call(rbind,lapply(rownames(tmat), function(idi){
    dfi <- data.frame(n=0,ntp=0)
    if(idi%in%rownames(u)){
      dfi$n <- sum(u[idi,])
      dfi$ntp <- sum(u[idi,]>0)
    }
    dfi
  }))
  
  
  x_opt2 <- cbind(data.frame(f_est=opt$par[1:nrow(tmat)],err=opt$value),dfstat)
  rownames(x_opt2) <- rownames(tmat)
  x_opt <- rbind(x_opt,x_opt2)
  return(x_opt)
}

##various functions for doing the polyharmonic splines:
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
  eE <- eigen(E)
  Uk <- eE$vectors[,1:k]
  Dk <- (diag(nrow(E))*eE$values)[1:k,1:k]
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
  k_range <- 5:min(100,length(y))
  #k_vals <- pbsapply(k_range,tst_k, x=x,f=y,nsplits=nsplits)
  k_vals <- sapply(k_range,function(ki){
    tryCatch(check_k(ki,x=x,f=y),error=function(e) return(Inf))
  })
  k_range[which.min(k_vals)]
}

optim_tps <- function(x_opt){
  x <- do.call(rbind,lapply(rownames(x_opt), function(r){
    as.numeric(unlist(strsplit(r,split="[.]")))
  }))
  y <- x_opt$f_est
  k <- find_best_k(x,y)
  ip <- eigenpoly(x,y,k)
  return(ip)
}

modify_config <- function(config,c.args){
  for(i in 1:nrow(c.args)){
    ix <- which(grepl(c.args$names[i],config))
    config[ix] <- paste(c.args$names[i],c.args$vals[i],sep=",")
  }
  return(config)
}

retrieve_fitness <- function(x_kary, x_times,x_path){
  ff <- paste0(x_path,stringr::str_pad(x_times,5,pad=0),".csv")
  fx <- lapply(ff,read.csv,header=F)
  nchrom <- ncol(fx[[1]])-2
  
  fitness <- unlist(sapply(fx, function(fxi){
    kary<-apply(fxi[,1:nchrom],1,paste,collapse=".")
    fitness<-fxi[,ncol(fxi)]
    names(fitness) <- kary
    fitness
  }))
  return(fitness[x_kary])
}
pij<-function(i, j, beta){
  qij <- 0
  if(abs(i-j)>i){ ## not enough copies for i->j
    return(qij)
  }
  # code fails for j = 0, but we can get this result by noting i->0 always
  ## accompanies i->2i
  if(j==0) j <- 2*i
  s <- seq(abs(i-j), i, by=2 )
  for(z in s){
    qij <- qij + choose(i,z) * beta^z*(1-beta)^(i-z) * 0.5^z*choose(z, (z+i-j)/2)
  }
  ## if there is a mis-segregation: ploidy conservation implies that both daughter cells will emerge simultaneously
  #if(i!=j) qij <- 2*qij
  
  return(qij)
}
getf <- function(pk,fobj){
  d <- apply(fobj$peaks,1,function(pki) sqrt(sum((pk-pki)^2)))
  sum(sin(d/fobj$wavelength)*fobj$scale)
}

gen_fitness_object <- function(cfig_path,lscape_path){
  lscape <- read.table(lscape_path,sep=",")
  cfig <- readLines(cfig_path)
  pk_scale <- cfig[grepl("scale",cfig)]
  pk_wl <- cfig[grepl("wavelength",cfig)]
  pk_scale <- as.numeric(tail(unlist(strsplit(pk_scale,split=",")),1))
  pk_wl <- as.numeric(tail(unlist(strsplit(pk_wl,split=",")),1))
  list(peaks=lscape,scale=pk_scale,wavelength=pk_wl)
}

gen_all_neighbours <- function(ids,as.strings=T){
  if(as.strings) ids <- lapply(ids, function(ii) as.numeric(unlist(strsplit(ii,split="[.]"))))
  nkern <- do.call(rbind,lapply(1:length(ids[[1]]), function(i){
    x0 <- rep(0,length(ids[[1]]))
    x1 <- x0
    
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0,x1)
  }))
  n <- do.call(rbind,lapply(ids, function(ii) t(apply(nkern,1,function(i) i+ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind,ids),n)
  n <- unique(n)
  n <- n[-(1:nids),]  
  n
}
