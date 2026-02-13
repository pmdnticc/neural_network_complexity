start_time <- Sys.time()
cat("Program started ... \n")
source("./cross_validation_v2.R")
source("./GDF_alt2.R")

library(nnet)

# function returns CV and GDF for every possible model for given data
# data.in: dataset to be analyzed
# var_names: variable names including predict and response variables
# nd: column # for the response variable in data.in
# kk: # of inversions in perturbation for responses
# nh: number of hidden units; dw: decay parameter in nnet settings
real_fit <- function(data.in,var_names,nd,kk=20,nh=2,dw=0.01,full=F){
  y_var <- var_names[1]
  xvars <- var_names[-1]
  nvar <- length(xvars)
  var_used <- pcv.m <- pcv.sd <- lcv.m <- lcv.sd <- gdf.m <- gdf.sd <- 
    setNames(replicate(nvar+1,array(),simplify=F), paste0("v",0:nvar))
  vlist <- var_used
  if(!full){
    pcv2.m <- pcv2.sd <- lcv2.m <- lcv2.sd <- lcv.m
  }
  for (i in 0:nvar) {
    vpos <- paste0("v",i)
    xx <- combn(xvars,i)
    var_used[[vpos]] <- xx
    for (j in 1:ncol(xx)) {
      form <- reformulate(xx[,j],response=y_var)
      model.fit <- do.call(nnet,list(
        formula=form,
        data=data.in,
        size=nh,
        linout=F,
        maxit=1000,
        entropy=T,
        decay=dw,
        trace=F
      ))
      if(full){
        tmp <- dlcv2(model.fit,data.in,nY=nd)
        pcv.m[[vpos]][j] <- tmp$ph
        pcv.sd[[vpos]][j] <- tmp$dph
        lcv.m[[vpos]][j] <- tmp$lcv
        lcv.sd[[vpos]][j] <- tmp$dlcv
        tmp <- get_GDF_h(model.fit,data.in,nY=nd,nk=kk,bootstrap=T)
        gdf.m[[vpos]][j] <- tmp$mean
        gdf.sd[[vpos]][j] <- tmp$sd
      } else{
        tmp <- dlcv2(model.fit,data.in,nY=nd)
        pcv2.m[[vpos]][j] <- tmp$ph
        pcv2.sd[[vpos]][j] <- tmp$dph
        lcv2.m[[vpos]][j] <- tmp$lcv
        lcv2.sd[[vpos]][j] <- tmp$dlcv
      }
    }
    if(!full){
      pos <- which.max(lcv2.m[[vpos]])
      pcv.m[[vpos]] <- pcv2.m[[vpos]][pos]
      pcv.sd[[vpos]] <- pcv2.sd[[vpos]][pos]
      lcv.m[[vpos]] <- lcv2.m[[vpos]][pos]
      lcv.sd[[vpos]] <- lcv2.sd[[vpos]][pos]
      form <- reformulate(var_used[[vpos]][,pos],response=y_var)
      model.fit <- do.call(nnet,list(
        formula=form,
        data=data.in,
        size=nh,
        linout=F,
        maxit=1000,
        entropy=T,
        decay=dw,
        trace=F
      ))
      tmp <- get_GDF_h(model.fit,data.in,nY=nd,nk=kk,bootstrap=T)
      gdf.m[[vpos]] <- tmp$mean
      gdf.sd[[vpos]] <- tmp$sd
      vlist[[vpos]] <- var_used[[vpos]][,pos]
    }
  }
  if(full){
    vlist <- var_used
  }
  list(lcv.m=lcv.m,lcv.sd=lcv.sd,pcv.m=pcv.m,pcv.sd=pcv.sd,
       gdf.m=gdf.m,gdf.sd=gdf.sd,var_list=vlist)
}



library(aplore3)
data(lowbwt)

set.seed(100)
# response bwt_low is binary (T,F)
lowbwt$bwt_low <- lowbwt$bwt<2500
nres <- dim(lowbwt)[2] # col number of the responses (T or F)
# all 8 prediction variables
vars_all <- c("bwt_low","age","lwt","race","smoke","ptl","ht","ui","ftv")

cat("Run function 'real_fit()' with different settings: \n")
cat("Hidden units nh = 2, 5, 10; and decay parameter dw = 0.01, 0.05, 0.1. \n")
for (hid in c(2,5,10)) {
  cat("Number of the hidden units: ", hid, "\n")
  for (dc in c(0.01,0.05,0.1)) {
    cat("Decay parameter: ", dc, "\n")
    start_internal <- Sys.time()
    tmp <- real_fit(data.in=lowbwt,var_names=vars_all,nd=nres,nh=hid,dw=dc)
    file_name <- paste0("real_fit_n",hid,"_dw",dc,"_new.rds")
    saveRDS(tmp,file=file.path("./data",file_name))
    end_internal <- Sys.time()
    print(end_internal-start_internal)
  }
}

end_time <- Sys.time()
cat("End of running the program \n")
print(end_time-start_time)

