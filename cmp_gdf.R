start_time <- Sys.time()
source("./GDF_alt2.R")
library(nnet)


set.seed(100)
n <- 200
N <- 100
px <- 0.5
py <- 0.3
pp <- 0.05
dw <- 0.1
hn <- 2
nx <- 8
xtype <- "mix" # options: binary, norm, mix
gdf.var <- gdf.vec <- array()
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
  data.rep <- data.frame(y,mget(paste0("x",1:nx)))
  model.fit <- nnet(form,data=data.rep,size=hn,linout=F,entropy=T,maxit=1000,
                    decay=dw,trace=F)
  
  # GDF horizontal
  yes <- rbinom(1,1,pp)==1
  if(yes){
    tmp <- get_GDF_h(model.fit,data.rep,nY=1,nk=20,bootstrap=T)
    gdf.vec[i] <- tmp$mean
    gdf.var[i] <- tmp$sd^2
  } else{
    gdf.vec[i] <- get_GDF_h(model.fit,data.rep,nY=1,nk=20)
    gdf.var[i] <- 0
  }
  
}

cat("GDF horizontal: ", mean(gdf.vec,na.rm=T), "+-", 
    sqrt(var(gdf.vec,na.rm=T)+mean(gdf.var,na.rm=T)),"\n")

end_time <- Sys.time()
print(end_time-start_time)


