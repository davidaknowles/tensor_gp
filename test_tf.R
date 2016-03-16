require(rstan)

source("load_data.R")

tenFactModel=stan_model("tensor_factorization_2.stan")
dat=list(K=3, D=length(levels(train$COMPOUND_A)), C=length(levels(train$CELL_LINE)), N=nrow(train), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), cellLine=as.numeric(train$CELL_LINE), synergy=train$SYNERGY_SCORE, Ntest=nrow(test), drugAtest=as.numeric(test$COMPOUND_A), drugBtest=as.numeric(test$COMPOUND_B), cellLinetest=as.numeric(test$CELL_LINE))
s=sampling(tenFactModel, data=dat, chains=7, cores=7, iter=1000, verbose=T)
samps=extract(s)

preds=colMeans(samps$test_synergy)
rmseOnTest(preds) # 33.7


