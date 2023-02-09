setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")
source("rscripts/polyh.R")
id <- "secondTest_rep_001"
config <- readLines(paste0("ABM/config/",id,".txt"))
dir <- "ABM/output/secondTest_rep_001/00000/proc_data/"

setwd(dir)
opt <- readRDS("optsimple.Rds")
opt <- opt[opt$err<1e6,]
plot(opt$f_est,opt$f_tru)

knots <- rownames(opt)
knots <- do.call(rbind,lapply(knots, function(ki){
  as.numeric(unlist(strsplit(ki,split="[.]")))
}))
f <- opt$f_est
k <- 2
ph <- setup_polyh(knots,f,k=k)
path <- "landscape.txt"

write_poly(path,knots,ph$w,ph$v)

config[grepl("fitness_landscape_type",config)] <- "fitness_landscape_type,polyh"
config[grepl("output_dir",config)] <- "output_dir,C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/proc_data/validation/"
config[grepl("fitness_landscape_file",config)] <- "fitness_landscape_file,C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/output/secondTest_rep_001/00000/proc_data/landscape.txt"
writeLines(config,"config.txt")

