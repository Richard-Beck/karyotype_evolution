
gaussf <- function(xy,pk1,pk2,sigma){
  max(exp(-sum((xy-pk1)^2/sigma^2)),
      exp(-sum((xy-pk2)^2/sigma^2)))
}

pk1 <- c(3,3)
pk2 <- c(7,7)
sigma <- 2

trn <- do.call(rbind,lapply(1:150, function(nothing){
  xy <- runif(2,min=0,max=10)
  f <- gaussf(xy,pk1,pk2,sigma)
  data.frame(x=xy[1],y=xy[2],f)
}))


library(caret)

set.seed(42)
# Set up the resampling, here repeated CV
tr <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# Note that trace is a parameter sent to the underlying modeling function
model <- train(f ~ ., data = trn, 
               #method = "lmStepAIC", 
               #method = "gaussprPoly",# good,limited
               #method = "krlsPoly",#  good here, not on actual data (needs more points??)
               #method = "rvmPoly",               
               #method = "svmPoly",
               #method="cubist",
               #method="gamboost",
               #method="treebag",
               trControl = tr)
model$results

tst <- expand.grid(x=seq(0,10,0.25),y=seq(0,10,0.25))
tst$f <- predict(model,tst)
tst$f_tru <- sapply(1:nrow(tst),function(i) gaussf(tst[i,1:2],pk1=pk1,pk2=pk2,sigma=sigma))

tst <- reshape2::melt(tst,id.vars=c("x","y"))

p <- ggplot(trn,aes(x=x,y=y,color=f))+
  geom_point()+
  scale_color_viridis_c()
p
p <- ggplot(tst,aes(x=x,y=y,fill=value))+
  facet_grid(cols=vars(variable))+
  geom_raster()+
  scale_fill_viridis_c()
p

