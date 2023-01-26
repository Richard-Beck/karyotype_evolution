library(deSolve)
library(parallel)



get_all_neighbours <- function(clone){
  
  do.call(rbind,lapply(1:length(clone), function(i){
    n_neighbours <- 
    do.call(rbind,lapply())
  }))
  
}


get_neighbour_fitness <-function(){
  clones <- names(x$x0)
  fitness <- x$fitness
  clones <- lapply(clones, function(ci) as.numeric(unlist(strsplit(ci,split="[.]"))))
  
  
  
}

runf <- function(fitness){
  
  mod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      m <- matrix(Pars,nrow=length(State))
      dx <- as.numeric(State%*%m)
      return(list(c(dx)))
    })
  }
  
  x0 <- x$x0
  names(x0) <- NULL
  mst <- x$d*fitness
  pars  <- c(mst=mst)
  yini <- c(x0=x0)
  times <- seq(0,ceiling(max(x$times)),1)
  out   <- ode(yini, times, mod, pars)
  out[,2:ncol(out)] <- out[,2:ncol(out)]/rowSums(out[,2:ncol(out)])
  out <- data.frame(out)
  colnames(out)[2:ncol(out)] <- rownames(x$x)
  
  
  
  return(out)
}


get_first_appearance <- function(dir,dt){
  ff <- list.files(dir)
  ff <- ff[!ff%in%c("log.txt","landscape.txt","proc_data")]
  ff <- ff[round(length(ff)*seq(0.1,0.8,0.1))]
  
  times <- ff
  times <- as.numeric(substr(times,1,5))
  times <- ceiling((times-times[1])*dt)
  i <- 1
  x <- do.call(rbind,lapply(1:length(ff), function(i) {
    fi <- ff[i]
    xi <- read.table(paste0(dir,"/",fi),sep=",")
    tt <- times[i]
    nchrom <- ncol(xi)-2
    fitness <- xi[,ncol(xi)]
    xi <- xi[1:nchrom]
    xi <- apply(xi,1,paste,collapse=".")
    tt <- rep(tt,length(xi))
    data.frame(karyotype=xi,gen=tt,fitness=fitness)
    }))
  
  x[!duplicated(x$karyotype),]
  
  
}

wrap_opt <- function(fitness) optim(par=fitness,fn=optimf,method = "L-BFGS-B",lower=0,upper=2)

read_landcape <- function(){
  landscape <- read.csv("landscape.txt",header=F)
  
  nchrom <- ncol(landscape)-2
  
  lapply(1:nrow(landscape), function(i){
    list( height =  landscape[i,nchrom+1],
          sigma  = landscape[i,nchrom+2],
          loc = landscape[i,1:nchrom])
  })
  
}

get_fitness <- function(clone,landscape){
  cx <- as.numeric(unlist(strsplit(clone,split="[.]")))
  max(sapply(landscape, function(li){
    li$height*exp(-sum(((li$loc-cx)/li$sigma)^2))
  }))
}

load_clones <- function(){
  ff <- list.files()
  ff <- ff[!ff%in%c("log.txt","landscape.txt")]
  ff <- ff[round(length(ff)*seq(0.1,0.8,0.1))]
  
  lapply(ff, function(fi) {
    xi <- read.csv(fi,header=F)
    nchrom <- ncol(xi)-2
    
    n <- xi[,nchrom+1]
    xi <- xi[,1:nchrom]
    xi <- apply(xi,1,paste,collapse=".")
    names(n) <- xi
    n
    
  })
}

load_fitness <- function(){
  ff <- list.files()
  ff <- ff[!ff%in%c("log.txt","landscape.txt")]
  ff <- ff[round(length(ff)*seq(0.1,0.8,0.1))]
  
  fi <- unlist(lapply(ff, function(fi) {
    xi <- read.csv(fi,header=F)
    nchrom <- ncol(xi)-2
    
    fi <- xi[,nchrom+2]
    xi <- xi[,1:nchrom]
    xi <- apply(xi,1,paste,collapse=".")
    names(fi) <- xi
    fi
    
  }))
  
  fi <- fi[!duplicated(names(fi))]
}

get_n <- function(yi,yt){
  
  n <- sapply(yt, function(yti){ 
    if(!yi%in%names(yti)) return(0)
    yti[yi]/sum(yti)
  })
  names(n) <- NULL
  return(n)
  
}

mod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    m <- matrix(Pars,nrow=length(State))
    dx <- as.numeric(State%*%m)
    return(list(c(dx)))
  })
}

optimf <- function(par,fit_dat){
  
  mod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      m <- matrix(Pars,nrow=length(State))
      dx <- as.numeric(State%*%m)
      return(list(c(dx)))
    })
  }
  fitness <- c(fitness,par)
  x0 <- fit_dat$pop[,1]
  names(x0) <- NULL
  mst <- fit_dat$d*fitness
  pars  <- c(mst=mst)
  yini <- c(x0=x0)
  meas_times <- as.numeric(colnames(fit_dat$pop))
  times <- seq(0,max(meas_times),1)
  out   <- ode(yini, times, mod, pars)
  out[,2:ncol(out)] <- out[,2:ncol(out)]/rowSums(out[,2:ncol(out)])
  meas <- out[out[,"time"]%in%meas_times,-1]
  meas <- t(meas)
  
  negll <- -log(as.numeric(unlist(lapply(1:nrow(meas), function(i) {
    dbinom(fit_dat$pop[i,],size=colSums(fit_dat$pop),prob=meas[i,])
  }
  ))))
  negll[!is.finite(negll)] <- 10^9
  
  negll <- sum(negll)
  print(negll)
  if(!is.finite(negll)){
    return(10^12)
  }
  return(negll)
  
}

est_clone_fitness <- function(zi,d,chr,fitness){
  
  n_neighbours <- 2*sum(as.numeric(unlist(strsplit(rownames(zi)[nrow(zi)],"[.]"))))
  medx <- sapply(rownames(zi)[-nrow(zi)], function(cli) med(rownames(zi)[nrow(zi)],cli,chr,string=T))
  medx <- (0.001/n_neighbours)^medx
  d <- rbind(d,0)
  d <- cbind(d,c(medx,1))
  
  fit_dat <- list(pop=zi,d=d,fitness=fitness)
  opt <- optimise(optimf,interval=c(0,3),fit_dat=fit_dat)
  opt$minimum
}



source("~/projects/008_birthrateLandscape/hoc_cpp/Rscripts/functions.R")
setwd("~/projects/008_birthrateLandscape/hoc_cpp/output/00000/")

log_file <- read.csv("log.txt")
landscape <- read_landcape()
x <- readRDS("proc_data/fit_dat.Rds")
clones <- readRDS("proc_data/clones.Rds")
cell_dat <- readRDS("proc_data/cell_dat.Rds")
z2 <- cell_dat[!rownames(cell_dat)%in%names(x$x0),]
z <- cell_dat[names(x$x0),]
dt <- log_file$dt

opt <- readRDS("proc_data/opt_out.Rds")
fitness <- opt$optim$bestmem

df <- runf(fitness)
df <- reshape2::melt(df,id.vars="time")
library(ggplot2)
y <- data.frame(t(x$x))
rownames(y) <- NULL
colnames(y) <- rownames(x$x)
y$time <- x$times
y <- reshape2::melt(y,id.vars="time")
p <- ggplot(df,aes(x=time,y=value))+
  facet_wrap(~variable)+
  geom_line()+
  geom_point(data=y)+
  theme(strip.text = element_text(size=6))+
  scale_y_sqrt()
p

#ggsave("proc_data/best.png",plot=p,width=10,height=7,units="in")

df_opt <- data.frame(estimated = opt$optim$bestmem,
                     true = x$fitness[names(x$x0)],
                     n = rowSums(z),
                     ntp = apply(z,1,function(zzi) sum(zzi>0)))

nz <- apply(z2,1,sum)
ntz <- apply(z2,1,function(zi) sum(zi>0))

z2 <- z2[nz>5&ntz>1,]

dz <- do.call(rbind,lapply(rownames(z2), function(x1){
  sapply(names(x$x0), function(x2) med(x1,x2,clones$chr,string = T))
}))

mindz <- apply(dz,1,min)

z2 <- z2[mindz==1,]
nz <- apply(z2,1,sum)
ntz <- apply(z2,1,function(zi) sum(zi>0))
z2 <- rownames(z2)
d0 <- x$d
fz2 <- sapply(z2,function(z2i){
  zi <- rbind(z,cell_dat[z2i,])
  rownames(zi)[nrow(zi)] <- z2i
  est_clone_fitness(zi,d0,clones$chr,fitness)
})
truz2 <- x$fitness[z2]

df2 <- data.frame(estimated=fz2,true=truz2,
                  n=nz,ntp=ntz)

df_opt$id <- "initial"
df2$id <- "secondary"

df_opt <- rbind(df_opt,df2)
df_opt$kary <- c(names(x$x0),z2)

df_plt <- df_opt[df_opt$n>5&df_opt$ntp>1,]

p <- ggplot(df_plt,aes(x=estimated,y=true,color=n,shape=id))+
  geom_point()+
  scale_color_viridis_c(trans="log")
p

saveRDS(df_opt,"proc_data/df_opt.Rds")
