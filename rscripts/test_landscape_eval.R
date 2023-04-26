  cellLine <- "SA609"
base.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(base.dir)
source("rscripts/analysis_functions.R")
source("rscripts/comparison_functions.R")
x0 <- readRDS(paste0("data/",cellLine,"/opt_data.Rds"))$x
x0 <- x0[apply(x0,1,max)>1,]
tt <- 1:ncol(x0)
sim.output.dir <- paste0("data/",cellLine,"/ABM_output/")
rna_data <- readRDS(paste0("data/",cellLine,"/opt_data.Rds"))

setwd(sim.output.dir)



ff <- list.files()

t0 <- as.numeric(sapply(ff, function(fi) tail(unlist(strsplit(fi,"_")),1)))
t1 <- t0+1

df <- do.call(rbind,lapply(1:7, function(i){
  print(i)
 
  ## get the cn distribution of the data at t=t0+1(passage)
  si <- list(x=do.call(rbind,lapply(rownames(x0), function(i) {
    as.numeric(unlist(strsplit(i,split="[.]")))
  })),
  n=x0[,t1[i]])
  si$x <- si$x[si$n>0,]
  si$n <- si$n[si$n>0]
  xi0 <- x0[,t1[i]]
  xit <- x0[,t0[i]]
  keep <- (xi0+xit)>0
  xi0 <- xi0[keep]
  xit <- xit[keep]
  r0ti <- xit/xi0
  
  
  
  dir <- paste0(ff[i],"/ABM_output/")
  ffi <- list.files(dir)
  do.call(rbind,lapply(1:length(ffi), function(counter) {
    diri <- paste0(dir,ffi[counter],sep="/") 
    
    ffij <- list.files(diri)
    ffij <- ffij[!ffij=="log.txt"]
    g <- as.numeric(substr(ffij,1,5)) ## available files in output dir
    
    dwass <- sapply(g, function(j){
      sj <- load_sim_dat(diri,j)
      get_dwass(si,sj)
    })
    
    errs <- sapply(g, function(j){
      xj <- proc_sim(diri,c(0,j))$x
      xitj <- xit[names(xit) %in% rownames(xj)]
      xj <- xj[names(xitj),]
      r0tj <- xj[,2]/xj[,1]
      r0tij <- r0ti[rownames(xj)]
      keep <- is.finite(r0tj) & is.finite(r0tij)
      mean((r0tj[keep]>1)==(r0tij[keep]>1))
      #R2(r0tj[keep],r0tij[keep])
    })
    
    
    Npass <- (1:length(dwass))
    
    gruup <- c(">1",">0",">3")[ceiling(counter/3)]
    
    data.frame(dwass,Npass,rep=counter,t0=t0[i],gruup,R2=errs)
  }))
  
}))
library(ggplot2)
p <- ggplot(df,aes(x=Npass,y=dwass))+
  facet_grid(cols=vars(t0),rows=vars(gruup),scales="free")+
  geom_line(aes(group=rep))+
  geom_smooth(aes(group=gruup,color=gruup))
p

p <- ggplot(df[df$gruup==">1",],aes(x=Npass,y=dwass))+
  facet_wrap(~t0,scales="free")+
  geom_line(aes(group=rep))+
  geom_smooth()+
  scale_x_continuous("arbitrary time unit")+
  scale_y_continuous("Wasserstein distance")
p


p <- ggplot(df,aes(x=Npass,y=R2))+
  facet_grid(cols=vars(t0),rows=vars(gruup),scales="free")+
  geom_line(aes(group=rep))+
  geom_smooth(aes(group=gruup,color=gruup))
p


