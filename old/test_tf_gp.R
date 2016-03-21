require(rstan)
require(doMC)
registerDoMC(7)

source("load_lb_data.R")

sm=stan_model("tensor_factorization_gp_3.stan")

cls=levels(train$CELL_LINE)

sqDist=dist$gex[cls,cls]^2

dat=list(K=3, D=length(levels(train$COMPOUND_A)), C=length(levels(train$CELL_LINE)), N=nrow(train), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), cellLine=as.numeric(train$CELL_LINE), synergy=train$SYNERGY_SCORE, Ntest=nrow(test), drugAtest=as.numeric(test$COMPOUND_A), drugBtest=as.numeric(test$COMPOUND_B), cellLinetest=as.numeric(test$CELL_LINE), sqDist=sqDist)
s=sampling(sm, data=dat, chains=7, cores=7, iter=1000, verbose=T)
samps=extract(s)

preds=colMeans(samps$test_synergy)
get_score(preds, test) # 0.0698
# 24.28925 RMSE

test$my_pred=preds
df=test %>% group_by(COMBINATION_ID) %>% summarise(pearson=cor(SYNERGY_SCORE, my_pred), n=length(SYNERGY_SCORE)) %>% as.data.frame

dfsub=df[df$n>2,]
ggplot(dfsub, aes(as.factor(n), pearson)) + geom_boxplot(outlier.shape = NA) + xlab("# cell lines in combination") + ylab("Pearson correlation") + geom_point(position = position_jitter(width = .5), size=3, alpha=.5) + theme_bw(base_size = 15)


