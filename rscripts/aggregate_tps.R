rbf <- function(r,eta=1e-9){
  if(r<eta) return(0)
  return(r^2*log(r))
}

## I think this does the rank reduction as explained in Wood (2003)
eigenpoly <- function(x,f,k){
  d <- as.matrix(dist(x))
  E <- do.call(rbind,lapply(1:nrow(d), function(i){
    sapply(d[i,],rbf)
  }))
  
  Tt <- cbind(1,x)
  
  
  ##potential issue - are we losing negative eigenvalues?
  ## R orders eigenvalues by size but perhaps it was supposed
  ## to be ordered by magnitude. 
  
  ##single value decomposition:
  eE <- eigen(E)
  Uk <- eE$vectors[,1:k]
  Dk <- (diag(nrow(E))*eE$values)[1:k,1:k]
  ##Ek:=Uk%*%Dk%*%Uk is the 'optimal' rank k approximation of E. 
  
  ##IDK what happens here
  A0 <- cbind(Uk%*%Dk,Tt)
  B <- t(Tt)%*%Uk
  z <- matrix(0,nrow(B),ncol(Tt))
  
  A <- rbind(A0,cbind(B,z))
  fe <- c(f,rep(0,nrow(B)))
  
  coef <- qr.solve(A,fe)
  
  ## wood says Uk%*%wk~=w where wk is the reduced coefficients.
  ## Therefore, presumably if we extract the approximate to w
  ## we can use this in the existing prediction routine. 
  
  list(coef=coef,w=Uk%*%coef[1:k],v=tail(coef,length(coef)-k),
       knots=x)
}

polypredict <- function(y,ip){
  r <- apply(ip$knots,1,function(i){
    di <- sqrt(sum((i-y)^2))
    rbf(di)
  })
  sum(c(r*ip$w,c(1,y)*ip$v))
}

check_k <- function(k,x,f){
  d <- as.matrix(dist(x))
  E <- do.call(rbind,lapply(1:nrow(d), function(i){
    sapply(d[i,],rbf)
  }))
  
  Tt <- cbind(1,x)
  
  
  ##potential issue - are we losing negative eigenvalues?
  ## R orders eigenvalues by size but perhaps it was supposed
  ## to be ordered by magnitude. 
  
  ##single value decomposition:
  Uk <- eigen(E)$vectors[,1:k]
  Dk <- (diag(nrow(E))*eigen(E)$values)[1:k,1:k]
  ##Ek:=Uk%*%Dk%*%Uk is the 'optimal' rank k approximation of E. 
  
  ##IDK what happens here
  A0 <- cbind(Uk%*%Dk,Tt)
  Ak <- A0%*%solve(t(A0)%*%A0)%*%t(A0)
  B <- t(Tt)%*%Uk
  z <- matrix(0,nrow(B),ncol(Tt))
  A <- rbind(A0,cbind(B,z))
  fe <- c(f,rep(0,nrow(B)))
  coef <- qr.solve(A,fe)
  fp <- A0%*%coef
  
  MSEk <- sum((f-fp)^2)
  tAk <- sum(1-diag(Ak))
  GCVk <- MSEk/(tAk/length(f))^2
  GCVk
}


find_best_k <- function(x,y){
  k_range <- 5:min(100,length(y))
  #k_vals <- pbsapply(k_range,tst_k, x=x,f=y,nsplits=nsplits)
  k_vals <- pbsapply(k_range,function(ki){
    tryCatch(check_k(ki,x=x,f=y),error=function(e) return(Inf))
  })
  k_range[which.min(k_vals)]
}

setup_validsims <- function(i){
  #create a sub directory for each train/test split
  sub.dir <- paste0(working.dir,stringr::str_pad(i,2,pad="0"),"/")
  dir.create(sub.dir)
  abm.out.dir <- paste0(sub.dir,"ABM_output")
  dir.create(abm.out.dir)
  setwd(paste0(root.dir,sim.dir))
  fo <- list.files()
  ##split initial sims into train and test
  f_train <- sample(fo,Ntrain)
  f_test <- fo[!fo%in%f_train]
  info <- list(f_train=f_train,f_test=f_test)
  
  saveRDS(info,paste0(sub.dir,"train_test.Rds"))
  ##write out the config file for the ABM to run
  init.config <- readLines(paste0(root.dir,"/config/",simname,".txt"))
  init.config[grepl("fitness_landscape_type",init.config)] <- "fitness_landscape_type,polyh"
  init.config[grepl("fitness_landscape_file",init.config)] <- paste0("fitness_landscape_file,",sub.dir,"landscape.txt")
  init.config[grepl("output_dir",init.config)] <- paste0("output_dir,",abm.out.dir)
  writeLines(init.config,paste0(sub.dir,"config.txt"))
  ##make batch file for running ABM
  cmd <- paste(cpp_source,paste0(sub.dir,"config.txt"))
  cmd <- sapply(1:Nreps,function(i) cmd)
  writeLines(cmd,paste0(sub.dir,"run_validation.bat"))
  
  ## estimate fitness landscape from training set
  df <- do.call(rbind,lapply(f_train, function(foi){
    readRDS(paste0(foi,"/proc_data/optsimple.Rds"))
  }))
  df <- df[df$err<1e3,]
  x <- do.call(rbind,lapply(rownames(df), function(r){
    as.numeric(unlist(strsplit(r,split="[.]")))
  }))
  y <- df$f_est
  k <- find_best_k(x,y)
  ip <- eigenpoly(x,y,k)
  saveRDS(ip,paste0(sub.dir,"polyscape.Rds"))
  
  ##write out the fitness landscape
  fscape <- cbind(ip$knots,ip$w)
  fscape <- rbind(fscape,ip$v)
  write.table(fscape,paste0(sub.dir,"landscape.txt"),row.names = FALSE,col.names=FALSE,sep=",")
  return(cmd)
}

library(pbapply)
Nreps <- 5 ## when running validation ABM sims, how many?
Nsplits <- 5 ## how many train/test splits
Ntrain <- c(1,3,5)
sims <- 1:5
cpp_source <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/bin/Debug/ABM.exe"
root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/"


for(sim in sims){
  cmds_sim <- c()
  for(nt in Ntrain){
    simname <- paste0("v00_rep_",stringr::str_pad(sim,3,pad=0))
    sim.dir <- paste0("output/",simname,"/")
    valid.dir <- paste0("validation/",simname,"/")
    out.dir <- paste0(root.dir,valid.dir)
    dir.create(out.dir)
    info <- paste0("Ntrain_",nt,"/")
    working.dir <- paste0(out.dir,info)
    dir.create(working.dir)
    cmds <- lapply(1:Nsplits, function(i){
      cmd <- tryCatch(setup_validsims(i),error=function(e) return(NULL))
      })
    cmds <- unlist(cmds)
    cmd_dir <- paste0(root.dir,"validation/cmds/",simname,".bat")
    cmds_sim <- c(cmds_sim,cmds) 
  }
  writeLines(cmds_sim,cmd_dir)
}













