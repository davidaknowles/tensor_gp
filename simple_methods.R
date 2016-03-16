source("load_data.R")

# Dumb method 0: predict 0
rmseOnTest(0) # 29.8

# Dumb method 1: predict global synergy mean.
sm=mean(train$SYNERGY_SCORE)
rmseOnTest(sm) # 26.5

# Method 2: predict average for drug combo
require(dplyr)
means=train %>% group_by(COMPOUND_A,COMPOUND_B) %>% summarize(mean=mean(SYNERGY_SCORE))
means=as.data.frame(means)
rownames(means)=paste(means$COMPOUND_A,means$COMPOUND_B,sep=".")
preds=means[test$COMBINATION_ID,"mean"]
rmseOnTest( preds ) # 24.1



