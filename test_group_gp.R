
require(rstan)
sm=stan_model("comb_therapy_models/group_gps.stan")

stopifnot(all( levels(mono$COMPOUND) == levels(train$COMPOUND_B) ))
stopifnot(all( levels(mono$COMPOUND) == levels(train$COMPOUND_A) ))
stopifnot(all( levels(mono$CELL_LINE) == levels(train$CELL_LINE) ))

cls=levels(train$CELL_LINE)
sqDist=lapply(dist,function(g) g[cls,cls]^2)

combinations=unique(train$COMBINATION_ID)

trainO=train[order(train$COMBINATION_ID),]
testO=test[order(test$COMBINATION_ID),]
s=table(train$COMBINATION_ID)
stest=table(test$COMBINATION_ID)
stopifnot(names(s) == names(stest))
dat=list(P=length(sqDist), N=nrow(trainO), Ntest=nrow(testO), C=length(cls), K=length(s), cellLines=as.numeric(trainO$CELL_LINE), cellLinesTest=as.numeric(testO$CELL_LINE), s=s, stest=stest, y=trainO$SYNERGY_SCORE, sqDist=sqDist )

o=optimizing(sm, dat, as_vector=F)

sqrt(mean( (testO$SYNERGY_SCORE-o$par$ytest)^2 ))

vb=rstan:::vb(sm, dat)
vm=getMeans(vb)
sqrt(mean( (testO$SYNERGY_SCORE-vm$ytest)^2 ))
