require(rstan)
source("load_data.R")

#stopifnot(all( levels(mono$COMPOUND) == levels(train$COMPOUND_B) ))
#stopifnot(all( levels(mono$COMPOUND) == levels(train$COMPOUND_A) ))
#stopifnot(all( levels(mono$CELL_LINE) == levels(train$CELL_LINE) ))

cls=levels(train$CELL_LINE)
sqDist=lapply(dist,function(g) g[cls,cls]^2)

all_dat=rbind(train,test[1:5,])

all_dat

dat=list(N=length(traininds), Ntest=length(testinds), y=y[traininds], C=length(y), P=2, sqDist=list( d, d2 ), cellLines=traininds, cellLinesTest=testinds )

o=optimizing(gpModel, data=dat,verbose=T,as_vector=F)

plot(x[traininds],y[traininds])
points( x[testinds], o$par$ytest, col="red", pch=16 )



