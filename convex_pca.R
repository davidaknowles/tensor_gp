# code to solve the optimization problem
# || X - Y ||_obs ^ 2 subject to ||X||_* < r
# where _obs denotes that we only care about the non-missing (NA) values in Y and
# ||X||_* is the nuclear (trace) norm of X

require(irlba)

convexPCA=function(Y,r,weights=NULL,Ytest=NULL,test_weights=NULL,its=1000,rmseTol=1e-3,warmStart=NULL,verbose=F){
  P=nrow(Y)
  N=ncol(Y)
  X=if (is.null(warmStart)) matrix(0,P,N) else warmStart
  if (is.null(weights)) weights=!is.na(Y)
  class(weights)="numeric"
  if (!is.null(Ytest)) {
    if (is.null(test_weights)) test_weights=!is.na(Ytest)
    class(test_weights)="numeric"
  }
  ertest=NA
  oldRmse=Inf
  for (it in 0:its){
    tol=max( 1e-1/(it+1)^2, 1e-6 )
    g=weights * (X-Y)
    svdTemp=irlba(-g,1,1,tol=tol)
    ruv=(r * svdTemp$u) %*% t(svdTemp$v)
    erv=(X-ruv) * weights
    stepSize=sum(g*erv)/sum(erv*erv)
    if (stepSize<0)
      cat('Warning: step size is',stepSize,'\n')
    stepSize=min(stepSize,1)
    X = (1.0-stepSize)*X + stepSize* ruv
    er=weights * (X-Y)
    if (!is.null(Ytest)) 
      ertest=(X-Ytest) * test_weights
    rmse=sqrt(sum(er^2)/sum(weights))
    rmseDelta=abs(rmse-oldRmse)
    if (verbose) 
      cat(it,rmse,sqrt(mean(ertest^2)),stepSize,rmseDelta,'\n')
    if ( rmseDelta < rmseTol)
      break
    oldRmse=rmse
  }
  X
}

cv.convexPCA=function(Y,weights=NULL,its=1000,rmseTol=1e-3) {
  
  if (is.null(weights)) weights=!is.na(Y)
  class(weights)="numeric"
  # break the matrix into training, test and validation sets, equally and at random
  rand=matrix( runif(nrow(Y)*ncol(Y)),nrow(Y))
  train=rand < 2/3
  test= !train & (weights>0)
  
  # training data
  train_weights=weights
  train_weights[!train]=0
  
  test_weights=weights
  test_weights[!test]=0
  
  # try a range of regularisation parameters
  rs=10^seq(0,log10(sum(Y^2,na.rm=T)+1),length.out = 10)
  trainErrors=numeric(length(rs))
  testErrors=numeric(length(rs))
  times=numeric(length(rs))
  for (ri in 1:length(rs)){
    r=rs[ri]
    cat("Trying regularisation parameter r=",r,"...\n")
    times[ri]=system.time( X<-convexPCA(Y,r,weights=train_weights,its=its,rmseTol=rmseTol,warmStart=if (ri>1) X else NULL ) )[1]
    trainErrors[ri]=sqrt(mean( train_weights * (X-Y)^2,na.rm=T))
    testErrors[ri]=sqrt(mean( test_weights * (X-Y)^2))
    cat("Error is",testErrors[ri],"\n")
  }
  rs[which.min(testErrors)]
}

wrapperPCA=function(Y,pcsRequired=5) {
  optimalR=cv.convexPCA(Y)
  cat("Picked regularisation parameter r=",optimalR,"\n")
  Yfull=convexPCA(Y,optimalR)
  irlba(Yfull,pcsRequired,pcsRequired)
}