require(data.table)

test=fread("../Drug_synergy_data/my_train_test_split/test.csv")
setDF(test)

train=fread("../Drug_synergy_data/my_train_test_split/train.csv")
setDF(train)

# TODO: check for overlaps between train and test
sum( duplicated( train[,c("CELL_LINE","COMPOUND_A","COMPOUND_B")] ) )

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