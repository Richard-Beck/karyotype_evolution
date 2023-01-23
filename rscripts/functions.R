## contains various functions for analysing karyotype data

##calculate the Mode of a vector x
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## wrapper for med function - returns a distance matrix
get_med <- function(y,chr){
  z <- do.call(rbind,lapply(strsplit(names(y),split="[.]"),as.numeric))
  
  do.call(rbind,lapply(1:nrow(z), function(i){
    sapply(1:nrow(z), function(j) med(z[i,],z[j,],chr))
  }))
}

## get the minimal spanning tree of a distance matrix
## 
get_mst <- function(d0,p0=0.01){
  dlog <- 1/p0^d0
  dtree <- ape::mst(dlog)
  d0*dtree
}

## function calculates the minimal event distance, i.e. least number
## of WGD's,depolyploidizations, whole chromosome or arm level mis-segregations
## needed to go from karyotype x1 to karyotype x2. 
med <- function(x1,x2,chr,string=FALSE){
  if(string){
    x1 <- as.numeric(unlist(strsplit(x1,split="[.]")))
    x2 <- as.numeric(unlist(strsplit(x2,split="[.]")))
  }
  ##Step 1. Keep proposing WGD's to cell x1 until it stops becoming more similar to x2
  d <- 0
  d0 <- sum(x1!=x2)
  while(sum((2*x1)!=x2)<d0){
    d <- d+1
    x1 <- 2*x1
    d0 <- sum(x1!=x2)
  }
  
  ##Step 1b keep proposing depolyploidizations 
  d0 <- sum(x1!=x2)
  while(sum((0.5*x1)!=x2)<d0){
    d <- d+1
    x1 <- 0.5*x1
    d0 <- sum(x1!=x2)
  }
  
  ##step 2. Do whole chromosome copy number changes to x1 until it stops becoming more similar to x2. 
  
  for(i in unique(chr)){
    i1 <- x1[chr==i]
    i2 <- x2[chr==i]
    
    while(mean(i1==i2)<1){
      ir <- pmin(i2/i1,2)
      d <- d+length(unique(ir[ir!=1]))
      i1 <- i1*ir
    }
  }
  
  return(d)
}

## this function gets all 1MS neighbours of a given karyotype. 
##KNOWN ISSUES:
##1) the function returns x1 as a neighbour of itself
##2) Unclear how to handle whole chromosome mis-segregations vs arm level
## missegregations. For e.g. assume chromosome 1p and 1q both have copy number
## of 2. Then it seems clear the possible 1MS neighbours are:
## (1,1); (3,3); (4;4) if whole chromosome MS, or (1,2); (3,2); (4,2);
## (2,1);(2,3);(2,4) if arm mis-segregation. Fine, but what are the 1 MS 
## neighbours if 1p has CN 1 and 1q has CN 3 (1,3) --> (?,?)
## the arm level 1MS neighbours are obvious but what does it mean 
## for this whole chromosome to mis-segregate once? I've decided
## a valid whole chrom mis-segregation is any allowable multiple
## of the initial state with all integer values, therefore
## (1,3) -> (2,6) only 

get_all_neighbours <- function(x1,chr){
  xn <- 2*x1
  xdpp <- 0.5*x1
  if(mean(xdpp==round(xdpp))==1) xn <- rbind(xn,xdpp)
  
  ##chrom arm neighbours
  xn <- rbind(xn,do.call(rbind,lapply(1:length(x1), function(i){
    do.call(rbind,lapply(1:(2*x1[i]), function(j){
      xij <- x1
      xij[i] <- j
      xij
    }))
  })))
  
  xn <- rbind(xn,do.call(rbind,lapply(unique(chr), function(i){
    i1 <- x1[chr==i]
    imin <- min(i1)
    do.call(rbind,lapply(1:(2*imin), function(j){
      xij <- x1
      xij[chr==i] <- (j/imin)*i1
      if(mean(xij==round(xij))!=1) return(NULL)
      return(xij)
    }))
  })))
  xn <- xn[!duplicated(xn),]
  return(xn)
}


## given two karyotypes x1,x2 return all the karyotypes that are 1MS
## away from x1 and are closer to x2 than x1 is. 
ms_step <- function(x1,x2,chr,string=F){
  if(string){
    x1 <- as.numeric(unlist(strsplit(x1,split="[.]")))
    x2 <- as.numeric(unlist(strsplit(x2,split="[.]")))
  }
  if(is.null(nrow(x1))) x1 <- matrix(x1,nrow=1)
  med0 <- min(apply(x1,1,function(i) med(i,x2,chr)))
  x1s <- do.call(rbind,lapply(1:nrow(x1), function(i) {
    get_all_neighbours(x1[i,],chr)
  }))
  is_closer <- apply(x1s,1,function(i) med(i,x2,chr)<med0) 
  x1s <- x1s[is_closer,]
  x1s[!duplicated(x1s),]
}