log <- read.csv("log.txt")
dt <- log$dt
ff <- list.files()
ff <- ff[!ff%in%c("log.txt","landscape.txt")]

output_dir <- "proc_data/"
dir.create(output_dir)




ff <- ff[round(length(ff)*seq(0.1,0.8,0.1))]

clones <- proc_clones(ff,Nclones = 3)
if(length(clones$yc)>1){
  
  ## there is some issue with the clones: 
  ## clones are appearing multiple times? 
  ## (there should be one entry per clone)  
  d <- get_med(clones$yc,clones$chr)
  dtree <- get_mst(d)

  
  if(max(dtree)>1){
    ##first only add the inbetween clones that are observed in the data
    ibtn <- gen_inbt(clones,dtree)
    clones$yc <- c(clones$yc,ibtn[ibtn%in%clones$y])
    clones$y <- clones$y[!clones$y%in%ibtn]
    
    ##check again if the tree is disconnected. If so, then we can add unobserved clones
    d <- get_med(clones$yc,clones$chr)
    dtree <- get_mst(d)
    if(max(dtree)>1){
      ibtn <- gen_inbt(clones,dtree)
      clones$yc <- c(clones$yc,ibtn)
      clones$y <- clones$y[!clones$y%in%ibtn]
    }
  }
  
  saveRDS(clones,paste0(output_dir,"clones.Rds"))
  
  ## measure the fitness of all observed clones
  x <- do.call(rbind,lapply(ff, read.table,sep=","))
  nchrom <- ncol(x)-2
  fitness <- x[,nchrom+2]
  x <- x[,1:nchrom]
  x <- apply(x,1,paste,collapse=".")
  names(fitness) <- x
  fitness <- fitness[!duplicated(names(fitness))]
  
  x <- lapply(ff, function(fi) {
    xi <- read.table(fi,sep=",")
    nchrom <- ncol(xi)-2
    n <- xi[,nchrom+1]
    xi <- xi[,1:nchrom]
    xi <- apply(xi,1,paste,collapse=".")
    names(n) <- xi
    yc <- clones$yc
    xi <- sapply(yc, function(ni) sum(n[names(n)==ni]))
    xii <- sapply(names(fitness), function(ni) sum(n[names(n)==ni]))
    
    list(xi=xi,xii=xii)
  })
  
  x2 <- do.call(cbind,lapply(x,function(xx) xx$xii))  
  x <-  do.call(cbind,lapply(x,function(xx) xx$xi))
  
  
  x0 <- x[,1]
  
  ## this part needs changed to be in line with GUSEV
  A <- do.call(rbind,lapply(clones$yc, function(yi){
    sapply(clones$yc, function(yj) wrap_pij(yi,yj,beta=0.001))
  }))
  
  times <- ff
  times <- as.numeric(substr(times,1,5))

  colnames(x2) <- times
  colnames(x) <- times
  
  fit_dat <- list(x0=x0,A=A,x=x,times=times,fitness=fitness)
  saveRDS(fit_dat,paste0(output_dir,"fit_dat.Rds"))
  saveRDS(x2,paste0(output_dir,"cell_dat.Rds"))
}



