library(transport)
library(ggplot2)

pop.comp <- function(g,f1,f2,info){
  x_test <- lapply(info$f_test, function(fti){
    fi <- paste0(f1,"/",fti)
    ffi <- list.files(fi)
    ffi <- ffi[!ffi%in%c("log.txt","proc_data")]
    xfi <- as.numeric(unlist(strsplit(ffi,split=".csv")))
    ffi <- ffi[which.min(abs(g-xfi))]
    ffi <- paste0(fi,"/",ffi)
    xi <- read.csv(ffi,header=F) 
    yi <- xi[,ncol(xi)-1]
    xi <- xi[,1:(ncol(xi)-2)]
    wpp(xi,yi/sum(yi))
  })
  
  x_pred <- lapply(list.files(paste0(f2,"ABM_output/")), function(fti){
    ffi <- list.files(fi)
    ffi <- ffi[!ffi%in%c("log.txt","proc_data")]
    xfi <- as.numeric(unlist(strsplit(ffi,split=".csv")))
    ffi <- ffi[which.min(abs(g-xfi))]
    ffi <- paste0(fi,"/",ffi)
    xi <- read.csv(ffi,header=F) 
    yi <- xi[,ncol(xi)-1]
    xi <- xi[,1:(ncol(xi)-2)]
    wpp(xi,yi/sum(yi))
  })
  
  ## calculate the mean wasserstein distance from x_test[i] to x_test[j!=i]
  ## and the mean dwass between x_test[i] and x_pred[j]
  
  
  df <- do.call(rbind,lapply(1:length(x_test), function(i){
    d_self <- mean(sapply((1:length(x_test))[-i], function(j){
      wasserstein(x_test[[i]],x_test[[j]])
    }))
    d_pred <- mean(sapply(x_pred, function(xj) wasserstein(x_test[[i]],xj)))
    data.frame(d_self,d_pred,g)
  }))
  return(df)
}

root.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/ABM/"
output.data.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/proc_dat/landscape_validation/"
output.fig.dir <- "~/projects/008_birthrateLandscape/karyotype_evolution/figures/"
sim.name <- "v00_rep_001"
output.data.dir <- paste0(output.data.dir,sim.name,"/")
dir.create(output.data.dir)
output.figure.name <- paste0(output.fig.dir,"wasserstein_",sim.name,".png")

g <- seq(200,2600,200)
Ntrain <- c(1,3,10)
Nreps <- 1:4
f1 <- paste0(root.dir,"output/",sim.name)
cmds <- expand.grid(g=g,Ntrain=Ntrain,Nreps=Nreps)

df <- do.call(rbind,pblapply(1:nrow(cmds), function(i){
  rep <- stringr::str_pad(cmds$Nreps[i],2,pad=0)
  Nt <- cmds$Ntrain[i]
  g <- cmds$g[i]
  f2 <- paste0(root.dir,"validation/",sim.name,"/Ntrain_",Nt,"/",rep,"/")
  info <- readRDS(paste0(f2,"train_test.Rds"))
  df <- do.call(rbind,lapply(g,pop.comp,f1=f1,f2=f2,info=info))
  df$rep<-rep
  df$Ntrain <- Nt
  df
}))


df$r <- df$d_pred/df$d_self

output.file.name <- paste0(output.data.dir,"wasserstein_ratio.Rds")
saveRDS(df,output.file.name)

df2 <- aggregate(list(ratio=df$r),by=list(g=df$g,
                                          rep=df$rep,
                                          Ntrain=df$Ntrain),
                 mean)

p <- ggplot(df2,aes(x=g,y=ratio))+
  facet_grid(rows=vars(Ntrain))+
  geom_point()+
  scale_x_continuous("timestep")+
  scale_y_continuous("Wasserstein ratio")
p



ggsave(output.figure.name,plot = p)





