source("~/projects/008_birthrateLandscape/karyotype_evolution/rscripts/regression_functions.R")
setwd("~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/")

predict_tps <- function(y,n,tpm,rot,cen){
  z <- ((y-cen)%*%rot)[1:n]
  print(z)
  dim(z) <- c(1,length(z))
  rich.predict2(z,predict.dat)
#  predict(tpm,z)
}

write_pd <- function(pd,dir,xpca){
  options(scipen=999)
  x <- c(paste0(c("nvar",ncol(pd$knots)),collapse=","),
         paste0(c("d",pd$d),collapse=","),
         paste0(c("xc",pd$xc),collapse=","),
         paste0(c("xs",pd$xs),collapse=","),
         paste0(c("c",pd$c),collapse=","),
         paste0(c("knots",c(t(pd$knots))),collapse=","),## laid out by row
         paste0(c("m",pd$null.args$m),collapse=","),
         paste0(c("p",c(pd$p)),collapse=","),
         paste0(c("cpca",c(xpca$center)),collapse=","),
         paste0(c("rpca",c(t(xpca$rotation))),collapse=",")) 
  writeLines(x,dir)
  
}

opt <- readRDS("proc_data/optsimple.Rds")
opt <- opt[opt$err<1e3,]
dx <- 1
x <- do.call(rbind,lapply(rownames(opt), function(oi) {
  xi <- as.numeric(unlist(strsplit(oi,split="[.]")))
  xi <- xi#*runif(length(xi),min=1-dx,max=1+dx)
  xi
}))
xraw <- x
xpca <- prcomp(x)

fits <- lapply(1:ncol(x), function(i){
  fit <- tryCatch({
    x <- xpca$x[,1:i]
    colnames(x) <- paste0("x",1:ncol(x))
    y <- opt$f_est
    Tps(x=x,Y=y)
  },error=function(e) return(NULL))
  return(fit)
})

n <- max(which(!sapply(fits,is.null)))

x <- xpca$x[,1:n]
colnames(x) <- paste0("x",1:ncol(x))
y <- opt$f_est
opt <- Tps(x=x,Y=y)
predict.dat <- list(null.function.name=opt$null.function.name,
                    null.args=opt$null.args,
                    xc=opt$transform$x.center,
                    xs=opt$transform$x.scale,
                    knots=opt$knots,
                    d=opt$d,c=opt$c,
                    ind.drift=opt$ind.drift,
                    p=opt$args$p)
founder <- rep(2,ncol(xraw))
predict_tps(founder,n,opt,xpca$rotation,xpca$center)
tst <- c(3,	1,	1,	1,	1,	2,	5,	1,	2,	1)
predict_tps(tst,n,opt,xpca$rotation,xpca$center)
write_pd(predict.dat,dir="proc_data/tst_scape.txt",xpca=xpca)