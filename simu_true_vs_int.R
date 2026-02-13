# helper: build true p given weights
expit <- function(z) exp(z)/(1+exp(z))

generate_data <- function(n=200, H=10, s=1, v0=-0.5) {
  # sample true weights (scale by s)
  vv <- rnorm(H,0,s)
  w0 <- rnorm(H,0,s)
  W <- matrix(rnorm(H*3,0,s),H,3)
  
  X <- matrix(rnorm(n*3),n,3)
  logits <- v0+rowSums(
    sapply(1:H,
           function(j) vv[j]*expit(
             w0[j]+W[j,1]*X[,1]+W[j,2]*X[,2]+W[j,3]*X[,3]
             )
           )
    )
  p <- expit(logits)
  y <- rbinom(n,1,p)
  list(X=X,y=y,p=p,logits=logits,
       diagnostics=c(mean_p = mean(p), var_p = var(p), mean_abs_logit = mean(abs(logits)),
                     bayes_acc = mean(pmax(p,1-p))))
}

# run for different scales
#scales <- c(0.5,1,2,3,5)
#res <- lapply(scales,function(s) generate_data(n=200,H=10,s=s,v0=1))
#sapply(res, function(x) x$diagnostics)

source("./cross_validation_v2.R")
source("./GDF_alt2.R")

library(nnet)
set.seed(100)
n_sample <- 200
pp <- 0.05 # we pick up a small subset of the total generated datasets for within data variance
N <- 100
h_unit <- 2
k1 <- 1
lcv.vec.var <- lcv.vec <- gdf.var <- gdf.vec <- array(dim=c(2,3,N))
lcv.mean <- lcv.sd <-gdf.sd <- gdf.mean <- array(dim=c(2,3))
lcv.paired <- gdf.paired <- array(dim=c(3,N))
k2 <- 1
for (dw in c(0.01,0.05,0.1)) {
  for (i in 1:N) {
    #if(i%%10==0) cat("# of replicates: ", i, "\n")
    data <- generate_data(n=n_sample,H=h_unit,s=1,v0=1)
    data_true.fit <- data.frame(y=data$y,x1=data$X[,1],x2=data$X[,2],x3=data$X[,3])
    y <- rbinom(n_sample,1,0.5)
    X <- matrix(rnorm(n_sample*3),n_sample,3)
    data_null.fit <- data.frame(y=y,x1=X[,1],x2=X[,2],x3=X[,3])
    
    #restore the random seed to ensure same initialization in nnet fit
    state <- .Random.seed 
    model_true.fit <- nnet(y~x1+x2+x3,
                           data=data_true.fit,size=h_unit,linout=F,
                           maxit=1000,entropy=T,decay=dw,trace=F,MaxNwts=1000)
    .Random.seed <- state
    model_null.fit <- nnet(y~x1+x2+x3, 
                           data=data_null.fit,size=h_unit,linout=F,
                           maxit=1000,entropy=T,decay=dw,trace=F,MaxNwts=1000)
    
    # same random seed for function dlcv2 under true and null
    state <- .Random.seed 
    tmp <- dlcv2(model_true.fit,data_true.fit,nY=1)
    lcv.vec[1,k2,i] <- tmp$ph
    lcv.vec.var[1,k2,i] <- tmp$dph^2
    .Random.seed <- state
    tmp <- dlcv2(model_null.fit,data_null.fit,nY=1)
    lcv.vec[2,k2,i] <- tmp$ph
    lcv.vec.var[2,k2,i] <- tmp$dph^2
    
    
    yes <- rbinom(1,1,pp)
      if(yes==1){
        state <- .Random.seed
        tmp <- get_GDF_h(model_true.fit,data_true.fit,nY=1,nk=20,bootstrap=T)
        gdf.vec[1,k2,i] <- tmp$mean
        gdf.var[1,k2,i] <- tmp$sd^2
        .Random.seed <- state
        tmp <- get_GDF_h(model_null.fit,data_null.fit,nY=1,nk=20,bootstrap=T)
        gdf.vec[2,k2,i] <- tmp$mean
        gdf.var[2,k2,i] <- tmp$sd^2
        } else {
        state <- .Random.seed
        gdf.vec[1,k2,i] <- get_GDF_h(model_true.fit,data_true.fit,nY=1,nk=20)
        .Random.seed <- state
        gdf.vec[2,k2,i] <- get_GDF_h(model_null.fit,data_null.fit,nY=1,nk=20)
        gdf.var[,k2,i] <- 0
        }
    }
  
  cat("Sample size: ", n_sample, "\n")
  cat("Hidden Units: ", h_unit, "\n")
  cat("Decay parameter: ", dw, "\n")
  
  cat("Under true model\n")
  cat("GDF: ", mean(gdf.vec[1,k2,]), "+-", 
      sqrt(var(gdf.vec[1,k2,])+mean(gdf.var[1,k2,])), "\n")
  cat("pCV: ", mean(lcv.vec[1,k2,]), "+-", 
      sqrt(var(lcv.vec[1,k2,])+mean(lcv.vec.var[1,k2,])), "\n")
  cat("Compare within and between variances: \n")
  cat("GDF within var: ", mean(gdf.var[1,k2,]), "between var: ", 
      var(gdf.vec[1,k2,]), "\n")
  cat("pCV within var: ", mean(lcv.vec.var[1,k2,]), "between var: ", 
      var(lcv.vec[1,k2,]), "\n")
  
  cat("Under intercept model\n")
  cat("GDF: ", mean(gdf.vec[2,k2,]), "+-", 
      sqrt(var(gdf.vec[2,k2,])+mean(gdf.var[2,k2,])), "\n")
  cat("pCV: ", mean(lcv.vec[2,k2,]), "+-", 
      sqrt(var(lcv.vec[2,k2,])+mean(lcv.vec.var[2,k2,])), "\n")
  cat("Compare within and between variances: \n")
  cat("GDF within var: ", mean(gdf.var[2,k2,]), "between var: ", 
      var(gdf.vec[2,k2,]), "\n")
  cat("pCV within var: ", mean(lcv.vec.var[2,k2,]), "between var: ", 
      var(lcv.vec[2,k2,]), "\n")
  
  # true model result
  lcv.mean[1,k2] <- mean(lcv.vec[1,k2,])
  lcv.sd[1,k2] <- sqrt(var(lcv.vec[1,k2,])+mean(lcv.vec.var[1,k2,]))
  gdf.mean[1,k2] <- mean(gdf.vec[1,k2,])
  gdf.sd[1,k2] <- sqrt(var(gdf.vec[1,k2,])+mean(gdf.var[1,k2,]))
  # null model result
  lcv.mean[2,k2] <- mean(lcv.vec[2,k2,])
  lcv.sd[2,k2] <- sqrt(var(lcv.vec[2,k2,])+mean(lcv.vec.var[2,k2,]))
  gdf.mean[2,k2] <- mean(gdf.vec[2,k2,])
  gdf.sd[2,k2] <- sqrt(var(gdf.vec[2,k2,])+mean(gdf.var[2,k2,]))
  
  # paired result
  lcv.paired[k2,] <- lcv.vec[1,k2,]-lcv.vec[2,k2,]
  gdf.paired[k2,] <- gdf.vec[1,k2,]-gdf.vec[2,k2,]
    
  k2 <- k2+1
}


sim_run <- list(
  lcv.mean=lcv.mean,
  lcv.sd=lcv.sd,
  gdf.mean=gdf.mean,
  gdf.sd=gdf.sd,
  lcv.paired=lcv.paired,
  gdf.paired=gdf.paired
)

saveRDS(sim_run,file="./data/run_n200_h2_x3.rds")





