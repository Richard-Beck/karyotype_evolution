Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

med <- function(x1,x2,chr,string=FALSE){
  if(string){
    x1 <- as.numeric(unlist(strsplit(x1,split="[.]")))
    x2 <- as.numeric(unlist(strsplit(x2,split="[.]")))
  }
  d <- 0
  ## was there a WGD?
  d1 <- sum(x1!=x2)
  dp <- Mode(x1/x2)
  x2i <- x2*dp
  
  
  d2 <- sum(x1!=x2i)+1
  ##if WGD, rescale x2 and add 1 to distance
  if(d2<d1){
    d <- d+1
    x2 <- x2i
  }
  
  ##how many other CNAs had to happen?
  d3 <- sapply(unique(chr), function(i){
    ## there are 4 cases
    ## 1) the entire chromosome is same for x1 and x2 (return 0)
    ## 2) there is 1 arm different between x1 and x2 (return 1)
    ## 3) whole chromosome is different (e.g. x1=2,2; x2=3,3) (return 1)
    ## 4) whole chromosome is different and 2 arm level CNAs must have (return 2)
    ## happened (e.g. x1=2,2; x2=1,3)
    i1 <- x1[chr==i]
    i2 <- x2[chr==i]
    
    ir <- unique(i1/i2) ## has length 1 if cases 1 & 3 are true, length 2 otherwise
    ir <- ir[!ir==1] ## removes one element if cases 2 or 3 are true
    length(unique(ir))
    
  })
  return(d+sum(d3))
}

## Define average number of cells with j chromosomes resulting from one cell with i chromosomes in the previous generation
pij<-function(i, j, beta){
  qij <- 0
  if(abs(i-j)>i){ ## not enough copies for i->j
    return(qij)
  }
  # code fails for j = 0, but we can get this result by noting i->0 always
  ## accompanies i->2i
  if(j==0) j <- 2*i
  s <- seq(abs(i-j), i, by=2 )
  for(z in s){
    qij <- qij + choose(i,z) * beta^z*(1-beta)^(i-z) * 0.5^z*choose(z, (z+i-j)/2)
  }
  ## if there is a mis-segregation: ploidy conservation implies that both daughter cells will emerge simultaneously
  #if(i!=j) qij <- 2*qij
  
  return(qij)
}

wrap_pij <- function(x1,x2,beta,string=TRUE){
  if(string){
    x1 <- as.numeric(unlist(strsplit(x1,split="[.]")))
    x2 <- as.numeric(unlist(strsplit(x2,split="[.]")))
  }
  
  prod(sapply(1:length(x1), function(i) pij(x1[i],x2[i],beta)))
  
}

proc_clones <- function(ff,Nclones=3){
  
  ## get the top N clones at each timepoint measured. How to choose N?
  yc <- unlist(lapply(ff, function(fi){
    xi <- read.table(fi,sep=",") 
    nchrom <- ncol(xi)-2
    yi <- xi[,1:nchrom]
    yi <- apply(yi,1,paste,collapse=".")
    n <- xi[,nchrom+1]
    names(n) <- yi
    n <- n[order(n,decreasing = T)]
    names(n)[1:Nclones]
  }))
  yc <- unique(yc)
  x <- do.call(rbind,lapply(ff, read.table,sep=","))
  nchrom <- ncol(x)-2
  
  n <- x[,nchrom+1]
  x <- x[,1:nchrom]
  y <- apply(x, 1, paste,collapse=".")
  
  y <- y[!y%in%yc]
  chr <- 1:ncol(x)
  list(yc=yc,y=y,chr=chr)
  
}

get_med <- function(y,chr,string=T){
  if(string) y <- do.call(rbind,lapply(strsplit(y,split="[.]"),as.numeric))

  do.call(rbind,lapply(1:nrow(y), function(i){
    sapply(1:nrow(y), function(j) med(y[i,],y[j,],chr))
  }))
}

mst_size <- function(d0){
  d <- 1/0.01^d0
  dtree <- ape::mst(d)
  sum(d*dtree)
}

get_mst <- function(d0,p0=0.01){
  dlog <- 1/p0^d0
  dtree <- ape::mst(dlog)
  d0*dtree
}

## given a tree and some clones, check for intermediate clones
check_inbetween <- function(dtree,clones){
  jumps <- which(lower.tri(dtree)&dtree>1)
  
  ids1 <- matrix(rep(names(clones$yc),length(clones$yc)), nrow=length(clones$yc),byrow = TRUE)
  ids2 <- matrix(rep(names(clones$yc),length(clones$yc)), nrow=length(clones$yc),byrow = FALSE)
  
  i <- ids1[jumps]
  j <- ids2[jumps]
  
  inbetweeners <- c()
  
  for(k in 1:length(i)){
    jump_size <- dtree[jumps[k]]
    
    d2 <- do.call(rbind,lapply(names(clones$y), function(y1){
      sapply(c(i[k],j[k]), function(y2) med(y1,y2,clones$chr,string=TRUE))
    }))
    rownames(d2) <- names(clones$y)
    
    worst_jumps <- apply(d2,1,max)
    inbetweeners <- c(inbetweeners ,rownames(d2)[worst_jumps<jump_size])
  }
  return(unique(inbetweeners))
}

## gets intermediates between karyotypes A and B assuming all chromosomes mis-segregate independently
## and no WGD. 
basic_inbt <- function(x1,x2,string=FALSE){
  if(string){
    x1 <- as.numeric(unlist(strsplit(x1,split="[.]")))
    x2 <- as.numeric(unlist(strsplit(x2,split="[.]")))
  }
  ms <- which(x1!=x2)
  ms_pos <- lapply(ms, function(i){
    xi <- rep(0,length(x1))
    xi[i] <- x2[i]-x1[i]
    return(xi)
  })
  
  
  
  do.call(rbind,lapply(1:(length(ms_pos)-1), function(i){
    co <- combn(1:length(ms_pos), m=i) 
    do.call(rbind,lapply(1:ncol(co), function(j){
      xi <- x1
      for(k in 1:nrow(co)){
        xi <- xi+ms_pos[[co[k,j]]]
      }
      xi
    }))
  }))
  
}

##generates all(*see basic inbt) in between karyotypes 
gen_inbt <- function(clones,dtree){
  jumps <- which(lower.tri(dtree)&dtree>1)
  ids1 <- matrix(rep(clones$yc,length(clones$yc)), nrow=length(clones$yc),byrow = TRUE)
  ids2 <- matrix(rep(clones$yc,length(clones$yc)), nrow=length(clones$yc),byrow = FALSE)
  
  i <- ids1[jumps]
  j <- ids2[jumps]
  
  ibtn <- unlist(lapply(1:length(i), function(k) {
    apply(basic_inbt(i[k],j[k],string=T),1,paste,collapse=".")
  }))
  unique(ibtn)
  
}
