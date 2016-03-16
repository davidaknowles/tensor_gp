
require(data.table)
require(Matrix)

ch1train=fread("../Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv")
setDF(ch1train)

ch1test=fread("../Drug_synergy_data/ch1_test_monoTherapy.csv")
setDF(ch1test)

ch1leader=fread("../Drug_synergy_data/ch1_leaderBoard_monoTherapy.csv")
setDF(ch1leader)

s=sort(table(c(ch1train$COMPOUND_A,ch1train$COMPOUND_B)),decreasing = T)[1:30]
par(mar=c(5,12,4,2)); barplot(s, horiz = T, las=2, xlab="# training pairs"); par(mar=c(5,4,4,2)); 

barplot(table(ch1train$QA), ylab="# training pairs", xlab="Quality metric")

allDrugs=union(ch1train$COMPOUND_A,ch1train$COMPOUND_B)

setdiff( union(ch1test$COMPOUND_A,ch1test$COMPOUND_B) ,allDrugs) # == NULL
setdiff(allDrugs, union(ch1test$COMPOUND_A,ch1test$COMPOUND_B)) # == NULL

setdiff( union(ch1leader$COMPOUND_A,ch1leader$COMPOUND_B) ,allDrugs) # == NULL

cellLines=unique(ch1train$CELL_LINE)

setdiff( cellLines, unique(ch1test$CELL_LINE) ) # "MDA-MB-175-VII" "HCT-116"       
setdiff( unique(ch1test$CELL_LINE), cellLines ) # NULL

setdiff( unique(ch1leader$CELL_LINE), cellLines ) # NULL

cellDrugFactor=function(ch1train) {
  ch1train$COMPOUND_A=factor(ch1train$COMPOUND_A, allDrugs)
  ch1train$COMPOUND_B=factor(ch1train$COMPOUND_B, allDrugs)
  ch1train$CELL_LINE=factor(ch1train$CELL_LINE,cellLines)
  ch1train
}

ch1train=cellDrugFactor(ch1train)
ch1test=cellDrugFactor(ch1test)
ch1leader=cellDrugFactor(ch1leader)

require(lattice)
pdf("train.pdf")
for (cl in cellLines){
  print(cl)
  submat=ch1train[ch1train$CELL_LINE==cl,]
  submat2=submat
  submat2$COMPOUND_A=submat$COMPOUND_B
  submat2$COMPOUND_B=submat$COMPOUND_A
  submat=rbind(submat,submat2)
  #sp=sparseMatrix( i=as.numeric(submat$COMPOUND_A), j=as.numeric(submat$COMPOUND_B), x=submat$SYNERGY_SCORE, dims=rep(length(allDrugs),2) ) 
  #image(sp)
  print(ggplot(submat, aes(COMPOUND_A,COMPOUND_B,fill=SYNERGY_SCORE)) + geom_tile() +  scale_fill_gradient2() + ggtitle(cl) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
}
dev.off()

require(dplyr)

plotPairs=function(ch1train) {
  submat2=ch1train
  submat2$COMPOUND_A=ch1train$COMPOUND_B
  submat2$COMPOUND_B=ch1train$COMPOUND_A
  submat=rbind(ch1train,submat2)
  counts=submat %>% group_by(COMPOUND_A,COMPOUND_B) %>% summarize(n=length(SYNERGY_SCORE))
  counts=as.data.frame(counts)
  ggplot(counts, aes(COMPOUND_A,COMPOUND_B,fill=n)) + geom_tile() +  scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ geom_abline(slope=1,intercept=1)
}

plotPairs(ch1train)
plotPairs(ch1test)
plotPairs(ch1leader)

ch1train$norm=sign(ch1train$SYNERGY_SCORE) * log(abs(ch1train$SYNERGY_SCORE))

qqnorm(ch1train$SYNERGY_SCORE,pch=16)
qqnorm(ch1train$SYNERGY_SCORE[ch1train$QA==1],pch=16)

ch1train$norm=sign(ch1train$SYNERGY_SCORE) * sqrt(abs(ch1train$SYNERGY_SCORE))
qqnorm(ch1train$norm[ch1train$QA==1],pch=16)

ch1trainHQ=ch1train[ch1train$QA==1,]

f=function(g) if (length(g)>5) abs(mean(g)/sd(g)) else NA
counts=ch1trainHQ %>% group_by(COMPOUND_A,COMPOUND_B) %>% summarize(n=length(SYNERGY_SCORE),mean=mean(SYNERGY_SCORE),sd=sd(SYNERGY_SCORE),cv=f(SYNERGY_SCORE))
counts=as.data.frame(counts)
hist(counts$cv,20)
counts=counts[counts$n>10,]
plot(counts$sd ~ counts$mean)

counts[ order(counts$sd,decreasing = T)[1:2], ]

syn=ch1train$SYNERGY_SCORE[ch1train$COMPOUND_A=="MAP2K_1" & ch1train$COMPOUND_B=="PIK3C"]
hist( syn, 20 )

ch1=list(train=ch1train, trainHQ=ch1trainHQ, test=ch1test, lead=ch1leader)

combs=lapply( ch1, function(g) paste(g$COMPOUND_A,g$COMPOUND_B,sep=" ") )

setdiff(combs$test,combs$train) # NULL set
setdiff(combs$lead,combs$train) # NULL set

ta=lapply(combs, table)
sharedCombs=Reduce(union, lapply(ta,names))
ta=lapply(ta, function(g) { temp=g[sharedCombs]; temp[is.na(temp)]=0; temp })
plot(ta$train+runif(length(sharedCombs)),ta$test+runif(length(sharedCombs)),pch=16,col=rgb(0,0,0,.3) , xlim=c(0,35), xlab="training count for pair", ylab="test count for pair")
points(ta$train+runif(length(sharedCombs)),ta$lead+runif(length(sharedCombs)), col=rgb(1,0,0,.6))

testIndices=unlist( foreach(comb=names(ta$trainHQ)) %do% {
  ss=strsplit(comb," ")[[1]]
  compA=ss[1]
  compB=ss[2]
  w=which(ch1trainHQ$COMPOUND_A==compA & ch1trainHQ$COMPOUND_B==compB)
  ntest=floor(length(w)/3)
  sample(w, ntest)
} )

trainIndices=setdiff(1:nrow(ch1trainHQ), testIndices)

myTrain=ch1trainHQ[trainIndices,]
myTest=ch1trainHQ[testIndices,]

#write.csv(myTest,file="../Drug_synergy_data/my_train_test_split/test.csv",quote=F, row.names = F)
#write.csv(myTrain,file="../Drug_synergy_data/my_train_test_split/train.csv",quote=F, row.names = F)

