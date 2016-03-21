x=seq(0,10,by=0.2)
y=sin(x) + rnorm(length(x))*.1

d=as.matrix(dist(x))^2

irrelevant_x=rnorm(length(x))
d2=as.matrix(dist(irrelevant_x))^2

traininds=sample(length(y),floor(length(y)/2))
testinds=setdiff(1:length(y),traininds)

dat=list(N=length(traininds), Ntest=length(testinds), y=y[traininds], C=length(y), P=2, sqDist=list( d, d2 ), cellLines=traininds, cellLinesTest=testinds )
gpModel=stan_model("comb_therapy_models/gp_mkl.stan")
o=optimizing(gpModel, data=dat,verbose=T,as_vector=F)

plot(x[traininds],y[traininds])
points( x[testinds], o$par$ytest, col="red", pch=16 )



