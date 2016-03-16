
x=seq(0,10,by=0.2)
y=sin(x) + rnorm(length(x))*.1

d=as.matrix(dist(x))^2

traininds=sample(length(y),floor(length(y)/2))
testinds=setdiff(1:length(y),traininds)

dat=list(N=length(traininds), Ntest=length(testinds), y=y[traininds], sqDist=d[traininds,traininds], sqDistTest=d[testinds,testinds], sqDistTrainTest=d[traininds,testinds] )
o=optimizing(gpModel, data=dat,verbose=T,as_vector=F)
cat("Lengthscale:",sqrt(o$par$inv_rho_sq),"\n")
s=sampling(gpModel, data=dat,  init=list(o$par), iter=100, verbose=T, chains=1)
samps=extract(s)
plot(sqrt(samps$inv_rho_sq), sqrt(samps$eta_sq), pch=16)
plot(sqrt(samps$eta_sq), sqrt(samps$sigma_sq), pch=16)
preds=colMeans(samps$testMean)
plot(x[testinds],preds,col="red"); points(x[traininds],y[traininds])
cor.test(y[testinds], preds)

# gp
gpModel=stan_model("gp_pred.stan")
phenos=c("IC50","H","Einf")
result=foreach (phe = phenos) %do% {
  pheA=paste0(phe,"_A")
  pheB=paste0(phe,"_B")
  x=c( as.character(train$CELL_LINE[train$COMPOUND_A=="AKT"]), as.character(train$CELL_LINE[train$COMPOUND_B=="AKT"]) )
  y=c( train[ train$COMPOUND_A=="AKT", pheA], train[train$COMPOUND_B=="AKT", pheB] )
  
  xtest=c( as.character(test$CELL_LINE[test$COMPOUND_A=="AKT"]), as.character(test$CELL_LINE[test$COMPOUND_B=="AKT"]) )
  ytest=c( test[test$COMPOUND_A=="AKT", pheA], test[test$COMPOUND_B=="AKT", pheB] )
  
  if (phe=="IC50") {
    y=scale(log(y))
    ytest=log(ytest)*attr(y,"scaled:scale") + attr(y,"scaled:center")
  }
  dat=list(N=length(y), Ntest=length(xtest), y=as.numeric(y), sqDist=gedist[x,x]^2, sqDistTest=gedist[xtest,xtest]^2, sqDistTrainTest=gedist[x,xtest]^2)
  
  o=optimizing(gpModel, data=dat,verbose=T,as_vector=F)
  cat("Lengthscale:",sqrt(o$par$inv_rho_sq),"\n")
  
  s=sampling(gpModel, data=dat, init=list(o$par), iter=100, verbose=T, chains=1)
  samps=extract(s)
  
  preds=colMeans(samps$testMean)
  cor(preds,ytest)
}

result=unlist(result)
names(result)=phenos


