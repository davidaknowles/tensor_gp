
#------------- LASSO --------------
require(glmnet)
monotest$lassoPred=NA
lassoRes=foreach (drug = levels(mono$COMPOUND)) %dopar% {
  monodrug=mono[mono$COMPOUND==drug,]
  g=cv.glmnet( t(ge[,monodrug$CELL_LINE]) , monodrug$Einf, parallel=F  )
  monotestdrug=monotest[monotest$COMPOUND==drug,]
  as.numeric( predict(g, t(ge[,monotestdrug$CELL_LINE])) )
}
names(lassoRes)=levels(mono$COMPOUND)
foreach(drug=levels(mono$COMPOUND)) %do% 
{ monotest[monotest$COMPOUND==drug,"lassoPred"] = lassoRes[[drug]] }

lassoRmse=rmseOnMonoTest(monotest$lassoPred)
#nn=stan_model("nn_real.stan")

plot(numFact, rmses, pch=16, xlab="number of factors", ylab="RMSE", type="b", ylim=c(33,36))
abline(rmseOnMonoTest(mean(mono$Einf)),0,col="red")
abline(lassoRmse,0,col="blue")

