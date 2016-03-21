source("load_data.R")

# Method 3: tensor factorization
require(rstan)

tenFactModel=stan_model("comb_therapy_models/tensor_factorization.stan")
dat=list(K=3, D=length(levels(train$COMPOUND_A)), C=length(levels(train$CELL_LINE)), N=nrow(train), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), cellLine=as.numeric(train$CELL_LINE), synergy=train$SYNERGY_SCORE, Ntest=nrow(test), drugAtest=as.numeric(test$COMPOUND_A), drugBtest=as.numeric(test$COMPOUND_B), cellLinetest=as.numeric(test$CELL_LINE))
s=sampling(tenFactModel, data=dat, chains=7, cores=7, iter=1000, verbose=T)
samps=extract(s)
preds=colMeans(samps$test_synergy)
rmseOnTest(preds) # 1000 iterations, 7 chains -> RMSE of 26.1

#"Tensor Fact"=26.1, "Tensor fact+GE"=24.3,
rmseVals=c("Predict 0"=29.8, "Predict mean"=26.5, "Mean of comb"=24.1, "Kernel regression"=23.85267, "Multitask KR"=22.22583, "With pathways"=21.98915)
ggplot(data.frame(x=factor(names(rmseVals),names(rmseVals)),y=rmseVals),aes(x,y))+geom_bar(stat="identity")+theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("method")+ylab("test set RMSE") 
#-------------- 
