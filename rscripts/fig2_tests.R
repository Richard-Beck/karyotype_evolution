setwd("~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/sweep_01/N_10_w_0p5_rep_00/")
source("~/projects/008_birthrateLandscape/karyotype_evolution/rscripts/analysis_functions.R")
df0 <- readRDS("~/projects/008_birthrateLandscape/karyotype_evolution/figure_data/fig2/sweep_1.Rds")
head(df0)
list.files()
df <- readRDS("optsimple_ntp_16.Rds")

optimx <- function(f,u,fit_dat){
  
  ux <- u/colSums(fit_dat$x)
  up <- 1-ux
  fp <- (fit_dat$pop.fitness-ux*f)/up
  
  tx <- as.numeric(colnames(fit_dat$x))
  for(i in (1+min(which(u>0))):length(ux)){
    t <- fit_dat$dt*diff(tx)[i-1]
    ux[i] <- ux[i-1]*exp(t*f)
    up[i] <- up[i-1]*exp(t*mean(fp[i],fp[i-1]))
    sx <- ux[i]+up[i]
    ux[i]<- ux[i]/sx
    up[i]<- up[i]/sx
  }
  negll <- -log(dbinom(u,colSums(fit_dat$x),prob=ux))
  negll[!is.finite(negll)] <- 10^9
  sum(negll)
}
optimsimple <- function(x){
  ntp <- apply(x$x,1,function(i) sum(i>0))
  to_fit <- which(ntp>1)
  x$x <- x$x[to_fit,]
  ix <- 1:nrow(x$x)
  
  opt <- do.call(rbind,lapply(ix, function(i){
    u <- x$x[i,]
    opt <- optimise(optimx,interval=c(0,1.2),u=u,fit_dat=x)
    data.frame(f_est=opt$minimum,err=opt$objective,
               n=sum(u),ntp=sum(u>0))
  }))
  rownames(opt) <- rownames(x$x)
  opt[opt$err<1e4,]
  
}



ntp <- 16

sim_path <- "train/00000/"
measure.times=round(seq(3000/ntp,3000,3000/ntp))
x <- proc_sim(sim_path,times=measure.times)


optim2 <- function(par,a,b){
  sum((par%*%a-b)^2)
}

a <- x$x

n <- apply(a,1,sum)
a <- a[n>1,]

a <- a/colSums(a)


b <- x$pop.fitness

p0 <- apply(a,1,which.max)
p0 <- b[p0]

opt <- optim(p0,optim2,a=a,b=b)

x_opt <- data.frame(f_est=opt$par,n=rowSums(a))
x_opt$f_tru <- retrieve_fitness(rownames(a), colnames(x$x), sim_path)

plot(x_opt$f_tru,x_opt$f_est)

x_opt <- optimsimple(x)
x_opt$f_tru <- retrieve_fitness(rownames(x_opt), colnames(x$x), sim_path)
x_opt$epc <- x_opt$err/x_opt$n
x_opt$etru <- abs(x_opt$f_est-x_opt$f_tru)
x_opt <- x_opt[order(x_opt$etru,decreasing=T),]

library(ggplot2)

p <- ggplot(x_opt,aes(x=f_tru,y=f_est,color=epc))+
  geom_point()+
  scale_color_viridis_c(trans="log")
p


