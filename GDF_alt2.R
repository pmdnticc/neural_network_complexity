# alternative algorithm for GDF version 2
# in this version the vertical method invert one component at one time

# the following two functions tackle the slope and std using bootstrap
# yp_pred: n x M matrix of fitted probabilities across M perturbation refits
# Yp:      n x M matrix of perturbed responses used in those refits
# both produced by your get_GDF_h() loop (before the lm step)
compute_gdf_slopes <- function(yp_pred,Yp,eps=1e-12) {
  n <- nrow(yp_pred)
  slopes <- numeric(n)
  for (i in 1:n) {
    # regression of predicted probs on perturbed responses for obs i
    # center both to improve numeric stability (optional)
    yhat <- yp_pred[i, ]
    ysim <- Yp[i, ]
    # small variance guard
    if (var(ysim)<eps) { 
      slopes[i] <- 0 
      } else {
        slopes[i] <- coef(lm(yhat ~ ysim))[2]
        }
  }
  sum(slopes)
}

bootstrap_gdf <- function(yp_pred,Yp,B=1000) {
  M <- ncol(yp_pred)
  gdf_boot <- numeric(B)
  for (b in 1:B) {
    idx <- sample.int(M,replace=T)  # resample columns
    gdf_boot[b] <- compute_gdf_slopes(yp_pred[,idx,drop=F],
                                      Yp[,idx,drop=F])
  }
  list(mean=mean(gdf_boot),
       sd=sd(gdf_boot))
}



# nk: number of perturbed flips; nn: number of replications
# nY: ordinal # of column containing responses from orig_data
get_GDF_h <- function(model_in,orig_data,nY,nk=1,nn=100,bootstrap=F,nB=100){
  # Horizontal Method
  Yo <- orig_data[,nY]
  nl <- length(Yo)
  Yp <- yp_pred <- data.frame(rep(0,nl))
  data <- orig_data
  k <- 1
  for (i in 1:nn) {
    # during each replication, nk perturbed flips are randomly selected
    flip <- rep(1,nl)
    while(length(which(flip==1))>0){
      to.flip <- which(flip==1)
      if(sum(flip)>nk){# ensure "to.flip" is an array not a scalar 
        flipped <- sample(to.flip,nk)
      }else{
        if(length(to.flip)>1) flipped <- sample(to.flip,sum(flip))
        else flipped <- to.flip
      }
      Yp[,k] <- Yo
      Yp[flipped,k] <- !Yo[flipped]
      data[,nY] <- Yp[,k]
      model_up <- update(model_in,data=data)
      yp_pred[,k] <- predict(model_up,type="raw")
      flip[flipped] <- 0 # label positions that are already performed perturbation
      k <- k+1
    }
  }
  yp_pred <- as.matrix(yp_pred)
  Yp <- as.matrix(Yp)
  if(bootstrap) {
    gdf <- bootstrap_gdf(yp_pred,Yp,B=nB)
  } else {
    gdf <- compute_gdf_slopes(yp_pred,Yp)
  }
  return(gdf)
#  hh.h <- lapply(1:length(Yo),function(j)lm(yp_pred[j,]~Yp[j,]))
#  gdf.h <- sum(sapply(hh.h,function(j)coefficients(j)[2]),na.rm=T)
#  return(gdf.h)
}




get_GDF_v <- function(model_in,orig_data,nY,nn=100){
  # Vertical Method
  Yo <- orig_data[,nY]
  nl <- length(Yo)
  hh.v <- array(dim=c(nl,nn))
  data <- orig_data
  for (i in 1:nn) {
    yo_pred <- predict(model_in,type="raw") # predictions have stochasticity
    # vertical method: only 1 component inverted for each
    flip <- rep(1,nl)
    while(sum(flip)>0){
      to.flip <- which(flip==1)
      flipped <- sample(to.flip,1)
      Ypp <- Yo
      Ypp[flipped] <- !Yo[flipped]
      data[,nY] <- Ypp
      model_up <- update(model_in,data=data)
      ypp_pred <- predict(model_up,type="raw")
      hh.v[flipped,i] <- 
        (ypp_pred[flipped]-yo_pred[flipped])/as.numeric(Ypp[flipped]-Yo[flipped])
      flip[flipped] <- 0
    }
  }
  list(
    mean=sum(hh.v)/nn,
    sd=sd(sapply(1:nn,function(k) sum(hh.v[,k])))
  )
#  gdf.v <- sum(hh.v)/nn
#  return(gdf.v)
}

