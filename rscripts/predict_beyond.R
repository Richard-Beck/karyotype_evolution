landscape <- read.table("landscape.txt",sep=",")
dat <- setup_opt("")
karys <- dat$trn[,-ncol(dat$trn)]

df <- make_nbrs(n=1:5,N=100,karys)
nn <- apply(df,1, function(dfi){
  min(as.matrix(dist(rbind(ii=dfi,karys)))[1,-1])
  
})
nn <- round(nn,digits=2)
f_nbrs <- apply(df,1,getf,landscape=landscape)

methods <- c("polyh","gaussprPoly","lmStepAIC","krlsPoly","ranger")

xi <- lapply(methods,fit_and_test,dat=dat,f_nbrs=f_nbrs,df=df,nn=nn)

metrics <- do.call(rbind,lapply(xi, function(xij) xij$metrics))
results <- do.call(rbind,lapply(xi, function(xij) xij$results))

saveRDS(metrics,"proc_data/landscape_metrics.Rds")
saveRDS(results,"proc_data/landscape_results.Rds")


