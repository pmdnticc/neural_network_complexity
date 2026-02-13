# Deviance function
DD <- function(y_obs,y_pred){
  ntot <- length(y_obs)
  n1 <- sum(y_obs)
  n0 <- ntot-n1
  lr <- sum(y_obs*log(y_pred)+(1-y_obs)*log(1-y_pred))
  lr <- lr+ntot*log(ntot)-n1*log(n1)-n0*log(n0)
  lr <- 2*lr
  return(lr)
}


library(nnet)


set.seed(100)
n <- 200
N <- 1000
px <- 0.5
py <- 0.3
dw <- 0.05
hn <- 10
nx <- 8
xtype <- "mix" # options: binary, norm, mix
dev <- array()
for (i in 1:N) {
  y <- rbinom(n,1,py)
  if(xtype=="binary") {
    if(nx==2){
      x1 <- rbinom(n,1,px)
      x2 <- rbinom(n,1,px)
    }
    if(nx==4){
      x1 <- rbinom(n,1,px)
      x2 <- rbinom(n,1,px)
      x3 <- rbinom(n,1,px)
      x4 <- rbinom(n,1,px)
    }
    if(nx==8){
      x1 <- rbinom(n,1,px)
      x2 <- rbinom(n,1,px)
      x3 <- rbinom(n,1,px)
      x4 <- rbinom(n,1,px)
      x5 <- rbinom(n,1,px)
      x6 <- rbinom(n,1,px)
      x7 <- rbinom(n,1,px)
      x8 <- rbinom(n,1,px)
    }
  }

  if(xtype=="norm") {
    if(nx==2){
      x1 <- rnorm(n)
      x2 <- rnorm(n)
    }
    if(nx==4){
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      x3 <- rnorm(n)
      x4 <- rnorm(n)
    }
    if(nx==8){
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      x3 <- rnorm(n)
      x4 <- rnorm(n)
      x5 <- rnorm(n)
      x6 <- rnorm(n)
      x7 <- rnorm(n)
      x8 <- rnorm(n)
    }
  }

  if(xtype=="mix") {
    if(nx==2){
      x1 <- rbinom(n,1,px)
      x2 <- rnorm(n)
    }
    if(nx==4){
      x1 <- rbinom(n,1,px)
      x2 <- rbinom(n,1,px)
      x3 <- rnorm(n)
      x4 <- rnorm(n)
    }
    if(nx==8){
      x1 <- rbinom(n,1,px)
      x2 <- rbinom(n,1,px)
      x3 <- rbinom(n,1,px)
      x4 <- rbinom(n,1,px)
      x5 <- rnorm(n)
      x6 <- rnorm(n)
      x7 <- rnorm(n)
      x8 <- rnorm(n)
    }
  }

  form <- reformulate(paste0("x",1:nx),response="y")
  model <- nnet(form,size=hn,linout=F,entropy=T,maxit=1000,decay=dw,trace=F)
  y_pred <- predict(model,type="raw")
  
  dev[i] <- DD(y,y_pred)
}


cat("mean LRT: ", mean(dev,na.rm=T), "\n")
