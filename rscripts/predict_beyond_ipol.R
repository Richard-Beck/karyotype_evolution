library(caret)

setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")

source("rscripts/regression_functions.R")
dir <- "ABM/output/secondTest_rep_003/00000/"
landscape <- read.table(paste0(dir,"landscape.txt"),sep=",")

make_fake_points <- function(kx){
  t(apply(kx,1,function(ki){
    i <- sample(1:length(ki),1)
    ki[i] <- 0
    ki
  }))
}
 
dat <- setup_opt(dir)


library(chebpol)


knots <- dat$trn[,1:(ncol(dat$trn)-1)]
f <- dat$trn$y
#f <- c(dat$trn$y,rep(0,nrow(knots)))
#knots <- rbind(knots,make_fake_points(knots))


model <- ipol(f,knots=apply(knots,1,as.numeric),method="poly")

m <- dat$tst[,1:(ncol(dat$tst)-1)]
f_est <- model(apply(m,1,as.numeric))

plot(f_est,dat$tst$y)

m <- dat$tst[,1:(ncol(dat$tst)-1)]
f_est <- model(apply(m,1,as.numeric))
plot(f_est,dat$tst$y)



results <- rbind(data.frame(true=dat$y_train_tru,estimated=dat$trn$y,id="a"),
                 data.frame(true=f_est,estimated=dat$tst$y,id="b"))


p <- ggplot(results,aes(x=estimated,y=true))+
  facet_grid(cols=vars(id))+
  geom_point()+
  geom_abline(color="red",size=1)
p


optimpeak <- function(par,landscape=NULL){
  if(!is.null(landscape)){
    return(-getf(par,landscape))
  }
  return(-model(par))
}

p0 <- rep(2,10)
optim(fn=optimpeak,par=p0,method = "L-BFGS-B",lower=0,upper=8)

x0 <- expand.grid(x1=seq(0,8,0.5),x2=seq(0,8,0.5))
f <- do.call(rbind,lapply(1:nrow(x0), function(i){
  xi <- rep(2,10)
  xi[c(3,7)] <- as.numeric(x0[i,])
  data.frame(f_tru=getf(xi,landscape = landscape),
             f_est = model(xi))
             
}))

x0 <- cbind(x0,f)
x0 <- reshape2::melt(x0,id.vars=c("x1","x2"))

p <- ggplot(x0,aes(x=x1,y=x2,fill=value))+
  facet_wrap(~variable)+
  scale_fill_viridis_c()+
  geom_raster()
p

wrap_ph <- function(pars){
  d1 <- 0
  d2 <- 0
  if(max(pars)>10) d1 <- max(pars)-10
  if(min(pars)<0) d2 <- 0-min(pars)
  return(-model(pars)+d1+d2)
} 


