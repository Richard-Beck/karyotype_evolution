
path1 <- "~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/"
path2 <- "~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/"

path2 <- "~/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/proc_data/validation/00000/"

f1 <- list.files(path1)
f1 <- f1[!f1%in%c("landscape.txt","log.txt","proc_data")]
f2 <- list.files(path2)
f2 <- f2[!f2%in%c("landscape.txt","log.txt","proc_data")]

t1 <- sapply(f1,function(fi) as.numeric(unlist(strsplit(fi,split=".csv"))))
t2 <- sapply(f2,function(fi) as.numeric(unlist(strsplit(fi,split=".csv"))))

tmax <- min(max(t1),max(t2))

t <- seq(0,tmax,500)

read.sim <- function(path) read.table(path,sep=",")

ft1 <- paste0(path1,sapply(t, function(ti) f1[which.min(abs(t1-ti))] ))
ft2 <- paste0(path2,sapply(t, function(ti) f2[which.min(abs(t2-ti))] ))


x1 <- lapply(ft1,read.sim)
x2 <- lapply(ft2,read.sim)

nchrom <- ncol(x1[[1]])-2



x <- rbind(do.call(rbind,x1),do.call(rbind,x2))
x <- x[,1:nchrom]
x <- unique(x)
tsne_ids <- apply(x,1,paste,collapse=".")

x1 <- do.call(rbind,lapply(1:length(x1),function(i){
  xi <- x1[[i]]
  kary <- apply(xi[,1:nchrom],1,function(cn) paste(cn,collapse="."))
  data.frame(kary,n=xi[,nchrom+1],f=xi[,nchrom+2],id="initial",t=t[i])
}))

x2 <- do.call(rbind,lapply(1:length(x2),function(i){
  xi <- x2[[i]]
  kary <- apply(xi[,1:nchrom],1,function(cn) paste(cn,collapse="."))
  data.frame(kary,n=xi[,nchrom+1],f=xi[,nchrom+2],id="validation",t=t[i])
}))

tsne <- data.frame(Rtsne::Rtsne(x)$Y)

x1 <- cbind(x1,tsne[sapply(x1$kary, function(id) which(tsne_ids==id)),])
x2 <- cbind(x2,tsne[sapply(x2$kary, function(id) which(tsne_ids==id)),])

y <- rbind(x1,x2)

library(ggplot2)
p <- ggplot(y,aes(x=X1,y=X2,color=n,size=n))+
  facet_grid(rows=vars(id),cols=vars(t))+
  geom_point()+
  scale_color_viridis_c()
p
