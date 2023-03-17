## downside of this augmentation is that we don't know the misseg rate.
## upside is that the model may not be particularly sensitive to this parameter.
## we may also be able to estimate the misseg rate in-pipeline.

root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/sweep_01/"
id <- "N_2_w_1p0_rep_02"
dir <- paste0(root.dir,id)  
setwd(dir)
sim_path <- "train/00000/"
cfig_path <- "config.txt"
lscape_path <- "landscape.txt"

fobj <- gen_fitness_object(cfig_path, lscape_path)

ntp <- 4
pm0 <- 0.0005

x <- proc_sim(sim_path,times=measure.times)

optimsimul <- function(par,y,dt,popf){
  y0 <- head(par,nrow(y))
  r0 <- abs(tail(par,nrow(y)))
  t <- as.numeric(colnames(y))
  lval <- do.call(rbind,lapply(1:length(y0), function(i){
    y0[i]+t*r0[i]*dt
  }))
  
  val <- apply(lval,2,function(lvi) {
    tmp <- exp(lvi-max(lvi))
    tmp/sum(tmp)
  })
  
  ll_rel <- -sum(unlist(lapply(1:nrow(y), function(i) {
    dbinom(y[i,],size=colSums(y), prob=val[i,],log=T)
  })))
  
  ll_vals <- -sum(dnorm(y0,mean=mean(popf),sd=sd(popf),log=T))
  
  ll_rel#+ll_vals
  
}

scalef <- function(s,estf,popf){
  sqrt(mean((s*estf-popf)^2))
}

ntp <- apply(x$x,1,function(i) sum(i>0))
n <- rowSums(x$x)
to_fit <- which(n>5)
y <- x$x[to_fit,,drop=FALSE]
y <- cbind(0,y)
colnames(y)[1] <- 0
nchrom <- ceiling(nchar(rownames(y)[1])/2)
founder_id <- paste(rep(2,nchrom),collapse=".")
if(!founder_id%in%rownames(y)){
  y <- rbind(0,y)
  rownames(y)[1] <- founder_id
}
y[founder_id,1] <- 1000
x_opt0 <- NULL
if(nrow(y)==1){
  x_opt0 <- data.frame(f_est=mean(x$pop.fitness),err=0,n=sum(y),ntp=ncol(y))
  rownames(x_opt0) <- rownames(y)
}else{
  dt <- x$dt
  # a good initial guess for fitness is pop. fitness when karyotype frequency peaks 
  
  
  
  r0 <- x$pop.fitness[apply(y[,-1],1,which.max)]
  t0 <- as.numeric(colnames(y)[apply(y,1,function(yi) which(yi>0)[1])])*dt
  y0 <- rnorm(nrow(y),sd=5)
  
  r0
  
  #y0 <- rep(-10,nrow(y))
  #nchrom <- ceiling(nchar(rownames(y)[1])/2)
  #founder_id <- paste(rep(2,nchrom),collapse=".")
  #y0[rownames(y)==founder_id] <- 0
  par <- c(y0,r0)
  popf<-x$pop.fitness
  opt <- optim(par,optimsimul,y=y,dt=dt,popf = popf)
  opts <- lapply(1:25, function(dumvar){
    y0 <- rnorm(nrow(y),sd=5)
    
    
    #y0 <- rep(-10,nrow(y))
    #nchrom <- ceiling(nchar(rownames(y)[1])/2)
    #founder_id <- paste(rep(2,nchrom),collapse=".")
    #y0[rownames(y)==founder_id] <- 0
    par <- c(y0,r0)
    
    opt <- optim(par,optimsimul,y=y,dt=dt,popf=x$pop.fitness)
  })
  
  best <- which.min(sapply(opts,function(oi) oi$value))
  opt <- opts[[best]]
  
  par <- abs(tail(opt$par,nrow(y)))
  names(par) <- rownames(y)
  
  estf <- par%*%apply(y[,-1],2,function(yi) yi/sum(yi))
  
  s <- optimise(scalef,interval=c(0,2),estf=estf,popf=x$pop.fitness)$minimum
  
  par <- par*s
  
  x_opt0 <- data.frame(f_est=par,err=opt$value,n=rowSums(y),ntp=apply(y,1,function(yi) sum(yi>0)))
  
}

pks <- lapply(rownames(x_opt0), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
x_opt0$f_tru <- sapply(pks,getf,fobj=fobj)


#v1 <- 0.001*exp((1:10))
#v2 <- exp((1:10)/3)
#lv1 <- log(0.001)+1:10
#lv2 <- log(1) +(1:10)/3
#v1/(v1+v2)
##since we are looking for R = exp(lv1)/(exp(lv1)+exp(lv2))
## we can rearrange: 1/R = (exp(lv1)+exp(lv2))/exp(lv1)
## 1/R = 1+exp(lv2-lv1)
#1/(1+exp(lv2-lv1))

#stop()


x_opt <- optimsimple(x)
pks <- lapply(rownames(x_opt), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
x_opt$f_tru <- sapply(pks,getf,fobj=fobj)
x_opt
x_opt0

x_opt$id <- "original"
x_opt0$id <- "updated"

xx <- rbind(x_opt,x_opt0)
library(ggplot2)
p <- ggplot(xx,aes(x=f_est,y=f_tru,color=id))+
  geom_point()+
  geom_abline()
p

stop()

xlim <- c(min(c(x_opt$f_est)))

plot(x_opt$f_est,x_opt$f_tru,xlim=c(-2,2))
lines(x_opt$f_tru,x_opt$f_tru)

plot(x_opt0$f_est,x_opt0$f_tru,xlim=c(-2,2))
lines(x_opt0$f_tru,x_opt0$f_tru)
print(opt$value)
stop()


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

optimf <- function(par,u,tmat,known_rates){
  
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

known_rates <- mean(x_opt$f_est)
par <- c(rep(0,nrow(tmat)),mean(x_opt$f_est),1)
opt <- optim(par,optimf,u=u,tmat=tmat,known_rates=known_rates)

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

pks <- lapply(rownames(x_opt), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
x_opt$f_tru <- sapply(pks,getf,fobj=fobj)
x_opt