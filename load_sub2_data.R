require(data.table)

train=fread("../Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv")
setDF(train)
train=train[train$QA==1,]

test=fread("../Drug_synergy_data/ch2_leaderBoard_monoTherapy.csv")
setDF(test)

allDrugs=Reduce(union, list(train$COMPOUND_A,train$COMPOUND_B,test$COMPOUND_A,test$COMPOUND_B))
cellLines=union(train$CELL_LINE,test$CELL_LINE)
cellDrugFactor=function(g) {
  g$COMPOUND_A=factor(g$COMPOUND_A, allDrugs)
  g$COMPOUND_B=factor(g$COMPOUND_B, allDrugs)
  g$CELL_LINE=factor(g$CELL_LINE,cellLines)
  g
}
train=cellDrugFactor(train)
test=cellDrugFactor(test)