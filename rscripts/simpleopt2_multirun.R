path <- "~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_002/"

get_pop_fitness <- function(foi){
  ff <- list.files(foi)
  ff <- ff[!ff%in%c("landscape.txt","log.txt","proc_data")]
  times <- ff
  times <- as.numeric(substr(times,1,5))
  dt <- read.csv(paste0(foi,"log.txt"),sep=",")$dt
  
  
  f <- sapply(paste0(foi,ff), function(fi){
    xi <- read.table(fi,sep=",")
    n <- xi[,ncol(xi)-1]
    f <- xi[,ncol(xi)]
    sum(n*f)/sum(n)
  })
  names(f) <- times
  f
}

getf <- function(t){
  t1 <- as.character(max(tf[as.numeric(names(f))<=t]))
  t2 <- as.character(min(tf[as.numeric(names(f))>=t]))
  mean(f[t1],f[t2])
}

optimx <- function(f,k,fit_dat){
  negll <- sapply(fit_dat, function(fdi){
    if(!k%in%row.names(fdi$xf)) return(0)
    ux <- fdi$xf[k,]/fdi$fit_dat$n
    up <- 1-ux
    fp <- (fdi$fit_dat$fx-ux*f)/up
    for(i in (1+min(which(ux>0))):length(ux)){
      t <- fdi$fit_dat$dt*diff(fdi$fit_dat$tx)[i-1]
      ux[i] <- ux[i-1]*exp(t*f)
      up[i] <- up[i-1]*exp(t*mean(fp[i],fp[i-1]))
      sx <- ux[i]+up[i]
      ux[i]<- ux[i]/sx
      up[i]<- up[i]/sx
      
    }
    negll <- -log(dbinom(fdi$xf[k,],fdi$fit_dat$n,prob=ux))
    negll[!is.finite(negll)] <- 10^9
    sum(negll)
  })
  sum(negll)
}

wrap_opt <- function(k,fit_dat){
  print(k)
  intervals <- list(c(1,2),c(0,3),c(0.5,2.8),c(0.5,3))
  ## apparently opt can miss if the function is not smooth enough and the interval is
  ## "unlucky". One way to avoid is by using a few differnt intervals - although maybe
  ## better in future if the Infinite values in optimx that cause this are not present
  opt <- lapply(intervals, function(i) optimise(optimx,interval=i,k=k,fit_dat=fit_dat))
  obj <- sapply(opt,function(oi) oi$objective)
  opt <- opt[[which.min(obj)]]
  data.frame(kary=k,f_est=opt$minimum,err=opt$objective)#,
  #n=sum(u),ntp=sum(u>0))
}

folders <- paste0(path,list.files(path)[1:10],"/")


opt_dat <- lapply(folders,function(foi) {
  x <- readRDS(paste0(foi,"proc_data/cell_dat.Rds"))
  tru_fitness <- readRDS(paste0(foi,"proc_data/fit_dat.Rds"))$fitness
  f <- get_pop_fitness(foi)
  list(x=x,f=f,tru_fitness = tru_fitness)
  })

tru_fitness <- unlist(sapply(opt_dat,function(oi) oi$tru_fitness))
tru_fitness <- tru_fitness[!duplicated(names(tru_fitness))]

fit_dat <- lapply(opt_dat,function(xi){
  x <- xi$x
  f <- xi$f
  Ncells <- colSums(x)
  ntp <- apply(x,1,function(i) sum(i>0))
  to_fit <- which(ntp>1)
  
  fx <- f[colnames(x)]
  tx <- as.numeric(names(fx))
  list(xf=x[to_fit,],fit_dat=list(n=Ncells,fx=fx,tx=tx,dt=dt,f=f),tru_fitness=xi)
})





unique_karyotypes <- unique(unlist(sapply(fit_dat, function(x) rownames(x$xf))))
df <- do.call(rbind,lapply(unique_karyotypes,wrap_opt,fit_dat=fit_dat))
df$f_tru <- tru_fitness[df$kary]
plot(df$f_tru,df$f_est)
