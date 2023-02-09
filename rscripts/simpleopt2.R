
getf <- function(t){
  t1 <- as.character(max(tf[as.numeric(names(f))<=t]))
  t2 <- as.character(min(tf[as.numeric(names(f))>=t]))
  mean(f[t1],f[t2])
}

optimx <- function(f,u,fit_dat){
  
  ux <- u/fit_dat$n
  up <- 1-ux
  fp <- (fit_dat$fx-ux*f)/up
  
  for(i in (1+min(which(u>0))):length(ux)){
    t <- fit_dat$dt*diff(fit_dat$tx)[i-1]
    ux[i] <- ux[i-1]*exp(t*f)
    up[i] <- up[i-1]*exp(t*mean(fp[i],fp[i-1]))
    sx <- ux[i]+up[i]
    ux[i]<- ux[i]/sx
    up[i]<- up[i]/sx
  }
  negll <- -log(dbinom(u,fit_dat$n,prob=ux))
  negll[!is.finite(negll)] <- 10^9
  sum(negll)
}

wrap_opt <- function(i,xf,fit_dat,fitness){
  u <- xf[i,]
  opt <- optimise(optimx,interval=c(0,1.2),u=u,fit_dat=fit_dat)
  
  data.frame(f_est=opt$minimum,f_tru=fitness[rownames(xf)[i]],err=opt$objective,
             n=sum(u),ntp=sum(u>0))
}


ff <- list.files()
ff <- ff[!ff%in%c("landscape.txt","log.txt","proc_data")]
times <- ff
times <- as.numeric(substr(times,1,5))
dt <- read.csv("log.txt",sep=",")$dt


f <- sapply(ff, function(fi){
  xi <- read.table(fi,sep=",")
  n <- xi[,ncol(xi)-1]
  f <- xi[,ncol(xi)]
  sum(n*f)/sum(n)
})
names(f) <- times


fitness <- readRDS("proc_data/fit_dat.Rds")$fitness
x <- readRDS("proc_data/cell_dat.Rds")
Ncells <- colSums(x)
ntp <- apply(x,1,function(i) sum(i>0))
to_fit <- which(ntp>1)

fx <- f[colnames(x)]
tx <- as.numeric(names(fx))
fit_dat <- list(n=Ncells,fx=fx,tx=tx,dt=dt,f=f)
xf<-x[to_fit,]
df <- do.call(rbind,lapply(1:nrow(xf),wrap_opt,xf=xf,fit_dat=fit_dat,fitness=fitness))
df$opt_err <- abs(df$f_est-df$f_tru)
rownames(df) <- rownames(xf)
saveRDS(df,"proc_data/optsimple.Rds")
