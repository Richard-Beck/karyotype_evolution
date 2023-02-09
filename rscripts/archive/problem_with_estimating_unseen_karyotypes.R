## situation 1

do_gen <- function(n,p,f=c(1,1,1)){
  n11 <- rbinom(1,n[1],f[1])
  n22 <- rbinom(1,n[2],f[2])
  n33 <- rbinom(1,n[3],f[3])
  n12 <- rbinom(1,n11,p)
  n13 <- rbinom(1,n11,p)
  
  n[1] <- n[1]+n11
  n[2] <- n[2]+n12+n22
  n[3] <- n[3] + n13+n33
  n
}

do_sim <- function(id,gens=250,p=0.01/80,f=c(0.9,0.9,1)){
  n <- c(10,0,0)
  n0 <- n
  for(i in 1:gens){
    if(sum(n)>10^7) n <- round(n/100)
    n <- do_gen(n,p,f)
    n0 <- rbind(n0,n)
  }
  n0 <- data.frame(n0)
  n0$gen <- 0:gens
  n0$f <- n0[,2]/(n0[,1]+n0[,2]+n0[,3])
  n0$id <- id
  return(n0)
}

x <- do.call(rbind,lapply(1:100,do_sim))

library(ggplot2)

p <- ggplot(x,aes(x=gen,y=f,group=id))+
  geom_line()+
  scale_x_continuous("generation")+
  scale_y_continuous("relative frequency of unseen neighbor")
p

## situation 2

do_gen <- function(n,p,f=c(1,1,1)){
  n11 <- rbinom(1,n[1],f[1])
  n22 <- rbinom(1,n[2],f[2])
  n33 <- rbinom(1,n[3],f[3])
  n12 <- rbinom(1,n11,p)
  n23 <- rbinom(1,n22,p)
  
  n[1] <- n[1]+n11
  n[2] <- n[2]+n12+n22
  n[3] <- n[3] + n33 + n23
  n
}

do_sim <- function(id,gens=150,p=0.01/40,f=c(0.7,0.8,1)){
  n <- c(10,0,0)
  n0 <- n
  for(i in 1:gens){
    if(sum(n)>10^7) n <- round(n/100)
    n <- do_gen(n,p,f)
    n0 <- rbind(n0,n)
  }
  n0 <- data.frame(n0)
  n0$gen <- 0:gens
  n0$f <- n0[,3]/(n0[,1]+n0[,3])
  n0$id <- id
  return(n0)
}

x1 <- do.call(rbind,lapply(1:100,do_sim,gens=350,f=c(0.7,0.7,1)))
x1$par <- "fit intermediate"


x2 <- do.call(rbind,lapply(1:100,do_sim,gens=350,f=c(0.7,0.3,1)))
x2$par <- "unfit intermediate"

x <- rbind(x1,x2)

p <- ggplot(x,aes(x=gen,y=f,group=interaction(par,id)))+
  geom_line(aes(color=par))+
  scale_x_continuous("generation")+
  scale_y_continuous("relative frequency of clone 3")+
  scale_color_discrete("")
p
 