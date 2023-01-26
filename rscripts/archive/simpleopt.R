
setwd("~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/")

ff <- list.files()
ff <- ff[!ff%in%c("landscape.txt","log.txt","proc_data")]
times <- ff
times <- as.numeric(substr(times,1,5))


f <- sapply(ff, function(fi){
  xi <- read.table(fi,sep=",")
  n <- xi[,ncol(xi)-1]
  f <- xi[,ncol(xi)]
  sum(n*f)/sum(n)
})

#f <- f[-1]

#plot(f[-1],log(100)/diff(times*dt))

#f <- log(100)/diff(times*dt) ## fitness estimate from growth rate
names(f) <- times


x <- readRDS("proc_data/fit_dat.Rds")
fitness <- x$fitness
x <- x$x
fx <- f[colnames(x)]
x <- x[1:min(nrow(x),length(fx)),]
x <- apply(x,2,function(xi) xi/sum(xi))
fitness <- fitness[rownames(x)]

optimf <- function(par){
  fx_est <- colSums(x*par)
  sqrt(mean((fx_est-fx)^2))
}

par <- rnorm(nrow(x),mean=1,sd=0.1)

opt <- optim(par,optimf,method = "L-BFGS-B",lower=0,upper=3)
par <- opt$par
cbind(par,fitness)
plot(par,fitness)
