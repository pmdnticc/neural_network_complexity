start_time <- Sys.time()
cat("Program started ... \n")
source("./cross_validation_v2.R")
source("./GDF_alt2.R")

library(nnet)
set.seed(100)
n <- 200
N <- 100  # Number of Monte Carlo replicates
px <- 0.5
py <- 0.3
pp <- 0.05
k_flip <- 20
nh <- 5
dw <- 0.01
pcv.var <- pcv.vec <- gdf.var.h <- gdf.vec.h <- array(dim=c(k_flip,N))
gdf.var.v <- gdf.vec.v <- array()
pcv.sd <- pcv.m <- gdf.sd.h <- gdf.m.h <- array()
for (i in 1:k_flip) {
  #cat("Number of inversions: ", i, "\n")
  for (j in 1:N) {
    if(j%%10==0) cat("Number of MC iterations: ", j, "\n")
    y <- rbinom(n,1,py)
    x1 <- rbinom(n,1,px)
    x2 <- rbinom(n,1,px)
    x3 <- rbinom(n,1,px)
    data.rep <- data.frame(y=y,x1=x1,x2=x2,x3=x3)
    model.fit <- nnet(y~x1+x2+x3,data=data.rep,size=nh,linout=F,maxit=1000,
                      entropy=T,decay=dw,trace=F)
    
    # p_CV mean and within-data variance
    #t0 <- Sys.time()
    tmp <- dlcv2(model.fit,data.rep,nY=1,kfold=i)
    pcv.vec[i,j] <- tmp$ph
    pcv.var[i,j] <- tmp$dph^2
    #t1 <- Sys.time()
    #cat("After cross-validation\n")
    #print(t1-t0)
    
    # GDF mean and within-data variance (using subset of replicated dataset)
    #t0 <- Sys.time()
    yes <- rbinom(1,1,pp)==1
    if(yes){
      tmp <- get_GDF_h(model.fit,data.rep,nY=1,nk=i,bootstrap=T)
      gdf.vec.h[i,j] <- tmp$mean
      gdf.var.h[i,j] <- tmp$sd^2
    } else{
      gdf.vec.h[i,j] <- get_GDF_h(model.fit,data.rep,nY=1,nk=i)
      gdf.var.h[i,j] <- 0
    }
    #t1 <- Sys.time()
    #cat("After GDF horrizontal\n")
    #print(t1-t0)
    if(i==1){ # compute GDF vertical only if one inversion at a time
      #t0 <- Sys.time()
      tmp <- get_GDF_v(model.fit,data.rep,nY=1)
      gdf.vec.v[j] <- tmp$mean
      gdf.var.v[j] <- tmp$sd^2
      #t1 <- Sys.time()
      #cat("After GDF vertical\n")
      #print(t1-t0)
    }
  }
  pcv.m[i] <- mean(pcv.vec[i,])
  pcv.sd[i] <- sqrt(var(pcv.vec[i,])+mean(pcv.var[i,]))
  gdf.m.h[i] <- mean(gdf.vec.h[i,])
  gdf.sd.h[i] <- sqrt(var(gdf.vec.h[i,]+mean(gdf.var.h[i,])))
}

gdf.m.v <- mean(gdf.vec.v)
gdf.sd.v <- sqrt(var(gdf.vec.v)+mean(gdf.var.v))

sim_run <- list(
  pcv.m = pcv.m,
  pcv.sd = pcv.sd,
  gdf.m.h = gdf.m.h,
  gdf.sd.h = gdf.sd.h,
  gdf.m.v = gdf.m.v,
  gdf.sd.v=gdf.sd.v
)

saveRDS(sim_run,file="./data/sim_null_S1.rds")

end_time <- Sys.time()
cat("End of running the program \n")
print(end_time-start_time)

