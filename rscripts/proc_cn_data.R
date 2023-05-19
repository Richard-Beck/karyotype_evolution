

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## x is the cn data file, arms is the arm loci file
call_cn <- function(x,arms,level="chrom"){
  ##proc cn file
  df <- do.call(rbind,lapply(rownames(x), function(xi) unlist(strsplit(xi,split="_"))))
  df <- data.frame(df[,1:2])
  colnames(df) <- c("chrom","chrompos")
  df$arm <- sapply(1:nrow(df), function(i){
    chrom <- df$chrom[i]
    chrompos <- df$chrompos[i]
    is_q <- arms[arms$chrom==chrom&arms$arm=="q","start"]<chrompos
    is_p <- !is_q
    
    c("q","p")[c(is_q,is_p)]  
    
  })
  
  x <- cbind(df[,c("chrom","arm")],x)
  xcolnames <- colnames(x)
  if(level=="chrom")  x <- split(x,f=x$chrom)
  if(level=="arm") x <- split(x,f=interaction(x$chrom,x$arm))
  lx <- sapply(x,nrow)
  x <- x[lx>0]
  
  x <- data.frame(do.call(rbind,lapply(x,  function(xi){
    y <- c(xi[1,1:2],as.numeric(apply(xi[,3:ncol(xi)],2,Mode)))
    y <- as.character(y)
    
  })))
  colnames(x) <- xcolnames
  
  return(x)
}

assign_to_tp <- function(x,meta){

  x <- x[!x$chrom%in%c("X","Y"),]
  x$chrom <- as.numeric(x$chrom)
  x <- x[order(x$chrom,x$arm),]
  nx <- colnames(x)[-c(1:2)]
  
  
  
  lookup=sapply(strsplit(nx,".", fixed = T),"[[",2)
  times <-  meta$timepoint[match(lookup,meta$library_id)]
  
  lookup <- data.frame(id=nx,times=times)
  
  lookup <- split(lookup,f=lookup$times)
  
  x <- lapply(lookup, function(i){
    xi <- x[,colnames(x)%in%i$id]
    rnxi <- colnames(xi)
    xi <- apply(xi,1,as.numeric)
    rownames(xi)<- rnxi
    return(xi)
  })
  
  times <- as.character(sapply(lookup, function(li) li$times[1]))
  names(x) <- times
  
  return(x)
}

##prepare data to fit sim
proc_exp <- function(x,fitness_file=NULL,dt=5){
  tps <- unlist(strsplit(names(x),split=".Rds"))
  tps <- as.numeric(substr(tps,start=2,stop=3))
  
  
  x <- lapply(x,function(xi){
    k <- apply(xi,1,paste,collapse=".")
    n <- 1
    dfi <- data.frame(karyotype=k,n=n)
    dfi <- aggregate(list(n=dfi$n),by=list(karyotype=dfi$karyotype),sum)
    dfi <- dfi[order(dfi$n,decreasing=T),]
  })
  
  
  unique.karyotypes <- unique(do.call(rbind,x)$karyotype)
  
  x <- do.call(rbind,lapply(as.character(unique.karyotypes), function(k){
    sapply(x, function(xi) sum(xi$n[xi$karyotype==k]))
  }))
  colnames(x) <- tps
  rownames(x) <- unique.karyotypes  
  rsx <- rowSums(x)
  x <- x[order(rsx,decreasing=T),]
  
  #fitness <- rep(NaN,ncol(x))
  #names(fitness) <- tps
  #if(!is.null(fitness_file)){
   # fitness_data <- readRDS(fitness_file)
    #fitness[substr(fitness_data$tp,start=2,stop=3)] <- fitness_data$x
  #}
  
  
  if(!is.null(fitness_file)){
    print("need to reinstate ability to include fitness file!")
  }else{
    fitness <- NULL
  }
  
  list(x=x,pop.fitness=fitness,dt=dt)
  
}

library(xlsx)
setwd("~/projects/008_birthrateLandscape/karyotype_evolution/salehi_data/")
##read inputs
arms <- readRDS("arm_loci.Rds")
meta=read.xlsx("41586_2021_3648_MOESM3_ESM.xlsx",sheetIndex = 1)

input_dir <- "01_raw_cn/"
output_dir <- "02b_arm_level/"
ff <- list.files(input_dir)
ff <- ff[grepl("SA609h",ff)]
for(fi in ff){
  
  raw_data_file <- paste0(input_dir,fi)
  output_filename <- head(unlist(strsplit(fi,split=".csv")),1)
  output_file <- paste0(output_dir,output_filename,".Rds")
  
  x <- read.csv(raw_data_file,row.names="bin_name")
  x <- call_cn(x,arms,level = "arm")
  x <- assign_to_tp(x,meta)
  x <- proc_exp(x)
  
  saveRDS(x,output_file)
}


