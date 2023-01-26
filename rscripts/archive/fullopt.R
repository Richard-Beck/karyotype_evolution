library(deSolve)
library(parallel)

optimf <- function(fitness,dt=0.1){
  
  mod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      m <- matrix(Pars,nrow=length(State))
      dx <- as.numeric(State%*%m)
      return(list(c(dx)))
    })
  }
  x0 <- fit_dat$pop[,1]
  names(x0) <- NULL
  mst <- fit_dat$d*fitness
  pars  <- c(mst=mst)
  yini <- c(x0=x0)
  meas_times <- as.numeric(colnames(fit_dat$pop))
  meas_times <- ceiling((meas_times-meas_times[1])*dt)
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

x <- readRDS("proc_data/fit_dat.Rds")
#y <- readRDS("proc_data/cell_dat.Rds")

#pop <- y[names(x$x0),]
#d <- x$d

fit_dat <- list(pop=x$x,d=x$A)

population_size <- 200
opt <- list()



library(DEoptim)
n <- length(x$x0)
init_pop <- do.call(rbind,lapply(1:population_size, function(dummy) rnorm(n,1,0.1)))

opt <- DEoptim(optimf,lower=rep(0,n),upper=rep(3,n),
               control = DEoptim.control(itermax=200,initialpop = init_pop,
                                         parallelType=1,packages="deSolve",parVar="fit_dat")
               )

#0.953243    1.074219    0.803691    1.084502    1.023107    0.441238    1.177920    1.185905    1.037048    1.131803    1.069948    1.210047    1.190214    1.135135    0.863024    1.265909    0.971423    1.251366    1.179070    1.080055    1.205432    0.893637    1.088385    1.108416    0.953226    1.057104    1.080691    1.210881    0.938151    0.994369    1.093663    1.090454    0.883043    1.069114    0.994474    1.127742    1.110912
#0.989918    1.013460    0.966507    1.025138    0.745490    0.674165    1.072386    1.102849    1.025912    1.138379    1.081897    1.024847    1.119317    1.039798    0.971431    1.117083    1.154093    0.922055    1.132700    1.027791    1.064194    0.866410    1.141735    0.937435    0.926148    1.138236    1.048527    0.892194    1.039277    0.959603    1.145289    0.984474    0.964165    0.929729    1.103170    1.085048    1.015183
saveRDS(opt,"proc_data/opt_out.Rds")

