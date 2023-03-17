root.dir <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/"
setwd(root.dir)
source("rscripts/functions.R")
cpp.out.dir <- "ABM/output/"
cpp_source <- "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/bin/Debug/ABM.exe"
config_path <- "ABM/config/"
ff <- list.files(cpp.out.dir)
batchname <- "randomTest"
simnames <- ff[grepl(batchname,ff)]

cmds <- lapply(simnames,function(si){
  cfile <- readLines(paste0(config_path,si,".txt"))
  cfile[grepl("fitness_landscape_type",cfile)] <- "fitness_landscape_type,polyh"
  p1 <- paste0(cpp.out.dir,si)
  reps <- list.files(p1)  
  cmd <- sapply(reps, function(ri){
    cfile_ri <- cfile
    fscape_path <- paste0(root.dir,cpp.out.dir,si,"/",ri,"/proc_data/estimated_landscape.txt")
    cfile_ri[grepl("fitness_landscape_file",cfile_ri)] <- paste0("fitness_landscape_file,",fscape_path)
    output_path <- paste0(root.dir,cpp.out.dir,si,"/",ri,"/proc_data/validation/")
    cfile_ri[grepl("output_dir",cfile_ri)] <- paste0("output_dir,",output_path)
    dir.create(output_path)
    cfile_path <- paste0(root.dir,cpp.out.dir,si,"/",ri,"/proc_data/config.txt")
    writeLines(cfile_ri,cfile_path)
    fitness_landscape_rpath <- paste0(root.dir,cpp.out.dir,si,"/",ri,"/proc_data/polyscape.Rds")
    cmd <- NULL
    cmd <- tryCatch({
      fitness_landscape <- readRDS(fitness_landscape_rpath)
      fscape <- cbind(fitness_landscape$knots,fitness_landscape$w)
      fscape <- rbind(fscape,fitness_landscape$v)
      write.table(fscape,fscape_path,row.names = FALSE,col.names=FALSE,sep=",")
      cmd <- paste(cpp_source,cfile_path)
    },error=function(e) return(NULL))
    return(cmd)
  })
})

cmds <- as.character(unlist(cmds))
cmds



