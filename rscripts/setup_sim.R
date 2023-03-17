options(scipen=999)
gen_lscape <- function(founder,Npeaks){
  peaks <- list()
  pk <- founder
  heights <- c()
  ht <- 1
  sigmas <- c()
  for(i in 1:Npeaks){
    ids <- sample(1:length(founder),2)
    pk[ids] <- sapply(pk[ids], function(j) sample(c(abs(0:(j-1)),(j+1):(2*j)),1))
    peaks <- c(peaks,list(pk))
    ht <- ht+runif(1,0.1,0.4)
    heights <- c(heights,ht)
    sigmas <- c(sigmas,runif(1,3,10))
  }
  
  sapply(1:length(peaks), function(i) paste(c(peaks[[i]],round(heights[i],digits=3),round(sigmas[i],digits=3)),collapse=","))
}




gen_randscape <- function(founder,Npeaks,scalef,wavelength=1){
  
  f0 <- 0
  ## if scalef <- 1/(pi*sqrt(30))
  ## then this would make the diploid cell in the top 10% fittest
  ## clones:
  while(f0<=0.4){
    pk <- lapply(1:Npeaks, function(i){
      pk <- sample(0:10,length(founder),replace=T)
    })
    
    d <- sapply(pk,function(ci) {
      sqrt(sum((founder-ci)^2))
    })
    f0=sum(sin(d/wavelength)*scalef)
    
  }
  print(f0)
  
  peak <- sapply(pk, paste,collapse="," )
  
}


gen_config <- function(output_dir="output",init_kary=rep(2,10),fitness_landscape_type="gaussian",
                       fitness_landscape_file="landscapes/landscape.txt",dt=0.1,
                       p=0.00005,Nsteps=3000,init_size=100000,scalef=1){
  
  c(paste(c("init_kary", init_kary),collapse=","),
    paste(c("fitness_landscape_type",fitness_landscape_type),collapse=","),
    paste(c("fitness_landscape_file", fitness_landscape_file),collapse=","),
    paste(c("dt",dt),collapse=","),
    paste(c("p", p),collapse=","),
    paste(c("Nsteps", Nsteps),collapse=","),
    paste(c("output_dir", output_dir),collapse=","),
    paste(c("init_size", init_size),collapse=","),
    paste(c("scale", scalef),collapse=","))
}

gen_replicates <- function(i,batchname,Nreps){
  sim_name <- paste0(batchname,"_rep_",stringr::str_pad(i,width=3,pad="0"))
  cpp_source <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/bin/Debug/ABM.exe"
  cpp_lscape_dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/landscapes/"
  cpp_config_dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/config/"
  
  fitness_landscape_file <- paste0(cpp_lscape_dir,sim_name,".txt")
  config_file <- paste0(cpp_config_dir,sim_name,".txt")
  
  cpp_output_dir <- paste0("C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/",sim_name,"/")
  dir.create(cpp_output_dir)
  
  
  
  founder <- rep(2,10)
  scalef <- 1/(pi*sqrt(30))
  lscape <- gen_randscape(founder,Npeaks=30,scalef)
  cfig <- gen_config(fitness_landscape_type = "random",
                     fitness_landscape_file = fitness_landscape_file,
                     output_dir = cpp_output_dir,scalef=scalef)
  writeLines(cfig,config_file)
  writeLines(lscape,fitness_landscape_file)
  
  cmd <- paste(cpp_source,config_file)
  rep(cmd,Nreps)
  
}

batchname <- "v00"
Nruns <- 5
Nreps <- 25

cmds <- unlist(lapply(1:Nruns,gen_replicates,batchname=batchname,Nreps=Nreps))
cmdpath <- paste0("C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/cmds/",batchname,".bat")
writeLines(cmds,cmdpath)



