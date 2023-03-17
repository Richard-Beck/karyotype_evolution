## downside of this augmentation is that we don't know the misseg rate.
## upside is that the model may not be particularly sensitive to this parameter.
## we may also be able to estimate the misseg rate in-pipeline.

root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/sweep_01/"
source("C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/rscripts/analysis_functions.R")
id <- "N_6_w_1p0_rep_07"
dir <- paste0(root.dir,id)  
setwd(dir)
sim_path <- "train/00000/"
cfig_path <- "config.txt"
lscape_path <- "landscape.txt"


fobj <- gen_fitness_object(cfig_path, lscape_path)
ntp <- 4
pm0 <- 0.00005
measure.times<-round(seq(0,3000,3000/(ntp-1)))
x <- proc_sim(sim_path,times=measure.times)
x_opt <- optimsimple(x)
x_opt <- infer_nn(x,x_opt,pm0)

pks <- lapply(rownames(x_opt), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
x_opt$f_tru <- sapply(pks,getf,fobj=fobj)
plot(x_opt$f_est,x_opt$f_tru)
x_opt