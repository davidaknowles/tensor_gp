require(rstan)

source("load_data.R")

sm=stan_model("tensor_factorization_gp_nn.stan")

cls=levels(train$CELL_LINE)

sqDist=dist$gex[cls,cls]^2

dat=list(K=3, K2=4, D=length(levels(train$COMPOUND_A)), C=length(levels(train$CELL_LINE)), N=nrow(train), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), cellLine=as.numeric(train$CELL_LINE), synergy=train$SYNERGY_SCORE, Ntest=nrow(test), drugAtest=as.numeric(test$COMPOUND_A), drugBtest=as.numeric(test$COMPOUND_B), cellLinetest=as.numeric(test$CELL_LINE), sqDist=sqDist)
s=sampling(sm, data=dat, chains=7, cores=7, iter=1000, verbose=T)

#v=vb(sm, data=dat)

samps=extract(s)

preds=colMeans(samps$test_synergy)
rmseOnTest(preds)
# 194.2439 haha