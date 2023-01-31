combine_data <- function(dirs){
  opt <- do.call(rbind,lapply(dirs, function(di){
    xi <- readRDS(paste0(di,"/proc_data/optsimple.Rds"))
    xi$karyotype<-rownames(xi)
    xi
  }))
  print(nrow(opt))
  opt <- opt[opt$err<10^6,]
  print(nrow(opt))
  opt <- aggregate(list(f_est=opt$f_est,f_tru=opt$f_tru,
                        n=opt$n,ntp=opt$ntp),
                   by=list(karyotype=opt$karyotype),mean)
  print(nrow(opt))
  opt
}

setup_opt <- function(dir){
  
  opt <- readRDS(paste0(dir,"proc_data/optsimple.Rds"))
  opt <- opt[opt$err<10^6,]
  f <- readRDS(paste0(dir,"proc_data/fit_dat.Rds"))$fitness
  
  x_train <- rownames(opt)
  
  x_test <- readRDS(paste0(dir,"proc_data/cell_dat.Rds"))
  x_test <- rownames(x_test)
  x_test <- x_test[!x_test%in%x_train]
  y_train <- opt$f_est
  y_train_tru <- opt$f_tru
  y_test <- as.numeric(f[x_test])
  
  trn <- data.frame(do.call(rbind,lapply(x_train, function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  })))
  trn$y <- y_train
  
  tst <- data.frame(do.call(rbind,lapply(x_test, function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  })))
  tst$y <- y_test
  
  list(tst=tst,trn=trn,y_train_tru=y_train_tru)
}

train_model <- function(trn,method = "gaussprPoly"){
  tr <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
  model <- train(y ~ ., data = trn, 
                 method=method,
                 trControl = tr)
}

make_nbr <- function(kary,n=1){
  i <- sample(1:length(kary),n,replace=T)
  for(ii in i)  kary[ii] <- kary[ii] + sample(c(-1,1),1)
  kary
}

make_nbrs <- function(n,N,karys){
  df <- data.frame(do.call(rbind,lapply(n, function(ni){
    i <- sample(1:nrow(karys),N,replace=T)
    do.call(rbind,lapply(i, function(ii) make_nbr(karys[ii,],ni)))
  })))

  k1 <- apply(karys,1,paste,collapse=".")
  k2 <- apply(df,1,paste,collapse=".")
  df <- df[!k2%in%k1,]
  df[!duplicated(df),]
}

get_fitness <- function(x1,x2,sigma,pk){
  pk*exp(-sum((x1-x2)^2/sigma^2))
}
getf <- function(k,landscape){
  max(sapply(1:nrow(landscape), function (i){
    k2 <- landscape[i,1:length(k)]
    sigma <- landscape[i,ncol(landscape)]
    pk <- landscape[i,ncol(landscape)-1]
    get_fitness(k,k2,sigma,pk)
  }))
}

fit_and_test <- function(method,dat,f_nbrs,df,nn){
  pred_b <- NULL
  pred_c <- NULL
  if(method=="polyh"){
    
    model <- ipol(dat$trn$y,
                  knots=apply(dat$trn[,-ncol(dat$trn)],1,as.numeric),
                  method="poly")
    m <- dat$tst[,1:(ncol(dat$tst)-1)]
    pred_b <- model(apply(m,1,as.numeric))
    pred_c <- model(apply(df,1,as.numeric))
  }else{
    model <- train_model(dat$trn,method=method)
    pred_b <- predict(model,dat$tst)
    pred_c <- predict(model,df)
  }
  
  
  dfa <- data.frame(predicted = dat$trn$y,
                    true = dat$y_train_tru,id="frequent clones")
  
  dfb <- data.frame(predicted=pred_b,
                    true=dat$tst$y,id="rare clones")
  
  dfc <- data.frame(predicted=pred_c,
                    true=f_nbrs,
                    id=paste("unobserved\nclones",nn))
  
  
  df <- rbind(dfa,dfb,dfc)
  df$method <- method
  
  metrics <- data.frame(rbind(postResample(dfa$predicted,dfa$true),
                              postResample(dfb$predicted,dfb$true),
                              postResample(dfc$predicted,dfc$true)))
  metrics$method <- method
  metrics$ids <- c("frequent clones","rare clones","unobserved clones")
  
  metrics <- reshape2::melt(metrics,id.vars=c("method","ids"))
  
  
  list(metrics=metrics,results=df)
}