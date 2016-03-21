require(data.table)
require(doMC)

load_data=function(train_fns, test_fn, dat_dir="Drug_synergy_data/") {
  # "Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv"
  # "Drug_synergy_data/ch1_leaderBoard_monoTherapy.csv"
  train=foreach(fn=train_fns, .combine = rbind) %do% {
    df=fread(paste0(dat_dir,fn))
    setDF(df)
    df
  }
  train=train[train$QA==1,]
  
  test=fread(paste0(dat_dir,test_fn))
  setDF(test)
  
  allDrugs=Reduce(union, list(train$COMPOUND_A,train$COMPOUND_B,test$COMPOUND_A,test$COMPOUND_B))
  cellLines=union(train$CELL_LINE,test$CELL_LINE)
  cellDrugFactor=function(g) {
    g$COMPOUND_A=factor(g$COMPOUND_A, allDrugs)
    g$COMPOUND_B=factor(g$COMPOUND_B, allDrugs)
    g$CELL_LINE=factor(g$CELL_LINE,cellLines)
    g
  }
  list( train=cellDrugFactor(train), test=cellDrugFactor(test) )
}
