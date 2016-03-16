require(data.table)

train=fread("../Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv")
setDF(train)

train_lb=fread("../Drug_synergy_data/ch1_LB.csv")
setDF(train_lb)


train=rbind(train,train_lb)

train=train[train$QA==1,]

test=fread("../Drug_synergy_data/ch1_test_monoTherapy.csv")
setDF(test)

allDrugs=union(train$COMPOUND_A,train$COMPOUND_B)
cellLines=union(train$CELL_LINE,test$CELL_LINE)
cellDrugFactor=function(g) {
  g$COMPOUND_A=factor(g$COMPOUND_A, allDrugs)
  g$COMPOUND_B=factor(g$COMPOUND_B, allDrugs)
  g$CELL_LINE=factor(g$CELL_LINE,cellLines)
  g
}
train=cellDrugFactor(train)
test=cellDrugFactor(test)