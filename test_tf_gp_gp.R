
require(rstan)
tf=stan_model("comb_therapy_models/tf_gp_gp.stan")

stopifnot(all( levels(mono$COMPOUND) == levels(train$COMPOUND_B) ))
stopifnot(all( levels(mono$COMPOUND) == levels(train$COMPOUND_A) ))
stopifnot(all( levels(mono$CELL_LINE) == levels(train$CELL_LINE) ))

cls=levels(mono$CELL_LINE)
dat=list(K=8, D=length(levels(train$COMPOUND_A)), C=length(levels(train$CELL_LINE)), N=nrow(train), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), cellLine=as.numeric(train$CELL_LINE), syn=train$SYNERGY_SCORE/100, Ntest=nrow(test), drugAtest=as.numeric(test$COMPOUND_A), drugBtest=as.numeric(test$COMPOUND_B), cellLinetest=as.numeric(test$CELL_LINE), Nmono=nrow(mono), drugMono=as.numeric(mono$COMPOUND), cellLineMono=as.numeric(mono$CELL_LINE), y=(mono$Einf-50)/10, sqDist=gedist[cls,cls]^2)
s=rstan:::vb(tf, data=dat, eta_adagrad=0.01, iter=1200, output_samples=30, seed=1)
m=getMeans(s)

rmseOnTest(m$testMean*100)

samps=extract(s)
preds=colMeans(samps$test_synergy)
rmseOnTest(preds) # 1000 iterations, 7 chains -> RMSE of 26.1

rmseVals=c("Predict 0"=29.8, "Predict mean"=26.5, "Mean of comb"=24.1, "Tensor Fact"=26.1, "Tensor fact+GE"=24.3)
ggplot(data.frame(x=factor(names(rmseVals),names(rmseVals)),y=rmseVals),aes(x,y))+geom_bar(stat="identity")+theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("method")+ylab("test set RMSE")
#-------------- 


