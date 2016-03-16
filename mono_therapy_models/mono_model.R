require(rstan)
source("load_data.R")


head(train)

getMono=function(train) {
  cols=c("COMPOUND","MAX_CONC","Einf")
  train[,c("CELL_LINE",paste0(cols,"_A"))]
  monoa=train[,c("CELL_LINE",paste0(cols,"_A"))]
  monob=train[,c("CELL_LINE",paste0(cols,"_B"))]
  colnames(monoa)=c("CELL_LINE",cols)
  colnames(monob)=c("CELL_LINE",cols)
  rbind(monoa,monob)
}
mono=getMono(train)
monotest=getMono(test)

rmseOnMonoTest=function(g) sqrt(mean( (monotest$Einf-g)^2 ))

rmseOnMonoTest(0) # 57.8
rmseOnMonoTest(mean(mono$Einf)) # 35.44



