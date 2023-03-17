root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
script_dir <- paste0(root.dir,"rscripts/")
source(paste0(script_dir,"comparison_functions.R"))
sweep_dir <- paste0(root.dir,"ABM/output/sweep_01")
setwd(sweep_dir)
f0 <- list.files()
g_range <- c(500,1000,2000)

R2 <- function(obs,pred){
  1-sum((pred-obs)^2)/sum((obs-mean(obs))^2)
}
RMSE <- function(obs,pred){
  round(sqrt(mean((obs-pred)^2)),digits=2)
}

df <- do.call(rbind,lapply(f0,function(id){
  print(id)
  df <- tryCatch({
    f1 <- paste0(sweep_dir,"/",id)
    
    id <- unlist(strsplit(id,"_"))
    ids <- id[c(2,4,6)]
    names(ids) <- id[c(1,3,5)]
    
    
    
    fi <- list.files(f1)
    landscapes <- fi[grepl("landscape_",fi)]
    cfigs <- fi[grepl("config_ntp_",fi)]
    
    cfigs_order <- unlist(strsplit(cfigs,".txt"))
    
    cfigs_order <- as.numeric(sapply(cfigs_order, function(cfi) tail(unlist(strsplit(cfi,"_")),1)))
    cfigs <- cfigs[order(cfigs_order)]
    
    if(length(cfigs)<1) return(NULL)
    s0 <- list(x=matrix(c(rep(2,ids["N"])),nrow=1,byrow=T),
               n=1)
    
    do.call(rbind,lapply(g_range,function(g){
      sj <- load_sim_dat(paste0(f1,"/train/00000/"),g)
      fk <- list.files(paste0(f1,"/train"))
      sk <- NaN
      if(length(fk)>1){
        sk <- sapply(fk[-1], function(fki){
          ski  <- load_sim_dat(paste0(f1,"/train/",fki,"/"),g)
          get_dwass(sj,ski)
        })
      }
      
      
      do.call(rbind,lapply(1:length(cfigs), function(i){
        cfig <- cfigs[i]
        ntp <- unlist(strsplit(cfig,".txt"))
        ntp <- tail(unlist(strsplit(ntp,"_")),1)
        f2 <- list.files(paste0(f1,"/test5"))[i]
        f2 <- paste0(f1,"/test5/",f2,"/")
        si <- load_sim_dat(f2,g)
        d_0_train <- get_dwass(s0,sj)
        d_0_test <- get_dwass(s0,si)
        d_train_test <- get_dwass(si,sj)
        
        
        optsimple.path <- paste0(f1,"/optsimple_ntp_",ntp,".Rds")
        optsimple <- readRDS(optsimple.path)
        optsimple <- optsimple[optsimple$err<1e3,]
        
        data.frame(Nchrom=ids["N"],wavelength=ids["w"],ntp=ntp,rep=tail(id,1),g=g,
                   d_train_train = mean(sk,na.rm=T),
                   p_err = 1-pnorm(log(d_train_test),mean(log(sk)),sd(log(sk))),
                   d_0_train,d_0_test,d_train_test,nknots=nrow(optsimple),
                   RMSE=RMSE(optsimple$f_tru,optsimple$f_est),
                   R2=R2(optsimple$f_tru,optsimple$f_est))
        
      }))
    }))
  },error=function(e) return(NULL))
  
  
  
  
  
}))


saveRDS(df,paste0(root.dir,"figure_data/fig2/sweep_1_v5.Rds"))




