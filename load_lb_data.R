require(data.table)

train=fread("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv")
setDF(train)

test=fread("Drug_synergy_data/ch1_LB.csv")
setDF(test)

# check for overlaps between train and test
# sum( duplicated( train[,c("CELL_LINE","COMPOUND_A","COMPOUND_B")] ) )

allDrugs=union(train$COMPOUND_A,train$COMPOUND_B)
cellLines=unique(train$CELL_LINE)
cellDrugFactor=function(g) {
  g$COMPOUND_A=factor(g$COMPOUND_A, allDrugs)
  g$COMPOUND_B=factor(g$COMPOUND_B, allDrugs)
  g$CELL_LINE=factor(g$CELL_LINE,cellLines)
  g
}
train=cellDrugFactor(train)
test=cellDrugFactor(test)

rmseOnTest=function(g) sqrt(mean( (test$SYNERGY_SCORE-g)^2 ))


#ge=read.csv("processed_data/gex_imputed_with_ccle.csv", check.names = F)
#rownames(ge)=ge[,1]
#ge[,1]=NULL
