# below is the code computing log-likelihood cross-validation

# k-fold cross-validation function: directly returns delta lcv asymptotically 
# equivalent to GDF
# default numbers for folds (kfold=10) and repeatation (Nrep=100)
dlcv2 <- function(model,orig_data,nY,kfold=10,Nrep=100){
  # the assumed data structure has nY-th column Y responses
  resY <- as.numeric(orig_data[,nY])
  # arg model is fitted model with orig_data; llm is calculated below
  y_pred <- predict(model,type="raw")
  llm <- sum(dbinom(resY,1,y_pred,log=T))
  Ns <- length(resY)
  # for binary data, prevalence in each fold should be kept approx equal
  ones <- which(resY==1)
  zeros <- which(resY==0)
  
  lik.rep <- numeric(Nrep)
  phat <- numeric(Nrep)
  for (i in 1:Nrep) {
    fold.pos <- array() # record position of i-th fold (i=1,...,k)
    fold.pos[ones] <- sample(rep(1:kfold,length.out=length(ones)))
    fold.pos[zeros] <- sample(rep(1:kfold,length.out=length(zeros)))
    # attach folds position to original dataset
    data.com <- cbind(orig_data,fold.pos)
    lik.fold <- numeric(kfold)
    for (j in 1:kfold) {
      model.train <- update(model,data=data.com[fold.pos!=j,])
      prob.pred <- predict(model.train,newdata=data.com[fold.pos==j,],type="raw")
      lik.fold.vec <- 
        dbinom(as.numeric(data.com[fold.pos==j,nY]),1,prob.pred,log=T)
      lik.fold.vec[which(is.infinite(lik.fold.vec))] <- NA
      lik.fold[j] <- sum(lik.fold.vec,na.rm=T)
    }
    lik.rep[i] <- sum(lik.fold)
    phat[i] <- llm-lik.rep[i]
  }
  list(ph=mean(phat,na.rm=T),dph=sd(phat,na.rm=T),
       lm=llm,lcv=mean(lik.rep,na.rm=T),
       dlcv=sd(lik.rep,na.rm=T))
}


