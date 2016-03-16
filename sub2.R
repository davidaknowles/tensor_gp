require(rstan)
require(doMC)
require(sna)
source("load_data.R")
source("load_sub2_data.R")

sm=stan_model("comb_therapy_models/gp_multitask_mkl.stan")

cls=levels(train$CELL_LINE)
sqDist=lapply(dist,function(g) g[cls,cls]^2)

drugs=allDrugs

pathways=read.csv("DREAM CHALLENGE TABLE DRUGS-TARGETS LOUKIA.csv",row.names = 1,check.names = F)
pathways=pathways[,colSums(pathways)>1]
ldrugs=do.call(rbind,strsplit(rownames(pathways),"-"))[,2]
sort(drugs[ !drugs %in% ldrugs])
sort( ldrugs[ ! ldrugs %in% drugs ] )

rownames(pathways)=ldrugs
pathways=pathways[drugs,]
pathways=as.matrix(pathways)

graph_dist=geodist(pathways,count.paths = F, inf.replace = 200)$gdist[1:nrow(pathways),1:nrow(pathways)]/2
graph_dist=graph_dist/median(graph_dist[upper.tri(graph_dist)])

co=read.table("pca_imputed_mono.txt",header=T)
mono_dist=as.matrix(dist(t(co)))
mono_dist=mono_dist/median(mono_dist[upper.tri(mono_dist)])
drugs[ ! drugs %in% rownames(mono_dist)] # a lot not done here

sqDist_dr=list( graph_dist^2 ) 

test_combs=unique(test$COMBINATION_ID)
topred=as.data.frame(matrix(NA, length(test_combs), length(cls)))
colnames(topred)=cls
topred$COMBINATION_ID=test_combs
require(reshape2)
m=melt(topred, id.vars="COMBINATION_ID")
test_compounds=do.call(rbind,strsplit(m$COMBINATION_ID,".",fixed = T))
m$COMPOUND_A=factor(test_compounds[,1],drugs)
m$COMPOUND_B=factor(test_compounds[,2],drugs)
colnames(m)[2]="CELL_LINE"

dat=list(N=nrow(train), Ntest=2, y=train$SYNERGY_SCORE, C=length(cls), D=length(levels(train$COMPOUND_A)), P=length(sqDist), sqDist_cl=sqDist, P_dr=length(sqDist_dr), sqDist_dr=sqDist_dr, cellLines=as.numeric(train$CELL_LINE), cellLinesTest=c(1,1), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), drugATest=c(1,1), drugBTest=c(1,1) )

o_init=optimizing(sm, data=dat,verbose=T,as_vector=F)

dat=list(N=nrow(train), Ntest=nrow(m), y=train$SYNERGY_SCORE, C=length(cls), D=length(levels(train$COMPOUND_A)), P=length(sqDist), sqDist_cl=sqDist, P_dr=length(sqDist_dr), sqDist_dr=sqDist_dr, cellLines=as.numeric(train$CELL_LINE), cellLinesTest=as.numeric(m$CELL_LINE), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), drugATest=as.numeric(m$COMPOUND_A), drugBTest=as.numeric(m$COMPOUND_B) )

o=optimizing(sm, data=dat,verbose=T, init=o_init$par, as_vector=F, iter=1)

names(o$par$eta_sq_cl)=names(sqDist)
barplot(sqrt(o$par$eta_sq_cl),ylab="importance")

#names(o$par$eta_sq_dr)=c("pathways","mono_therapy")
#barplot(sqrt(o$par$eta_sq_dr),ylab="importance")

N=nrow(train)
o$par$sigma_sq
SN=o$par$Sigma_no_noise

ytrain = SN %*% ( solve(SN + diag(N) * o$par$sigma_sq, dat$y - o$par$mu) ) + o$par$mu
train$err= (ytrain-dat$y)^2
plot(dat$y, ytrain)

require(dplyr)
errs=train %>% group_by(COMBINATION_ID) %>% summarize(err=mean(err), mean=mean(SYNERGY_SCORE))
setDF(errs)
a=as.data.frame(errs[1:167,])

plot( a$mean, sqrt(a$err ), xlab="mean synergy", ylab="training RMSE", pch=16 ) 
#text( a$mean, sqrt(a$err ), a$COMBINATION_ID)
a$cv=sqrt(a$err)/ abs(a$mean) 
a=a[order(a$cv),]
cv=a$cv
names(cv)=a$COMBINATION_ID
cv=cv[c(1:20,length(cv)-20+(1:20))]

b=a[c(1:20,nrow(a)-20+(1:20)),]
b$COMBINATION_ID=factor(b$COMBINATION_ID,b$COMBINATION_ID)
ggplot(b, aes(x=COMBINATION_ID, y=cv) ) + geom_bar(stat = "identity")+coord_flip() + scale_y_log10()

str(o)

rownames(a)=a$COMBINATION_ID
a[as.character(sub2$COMBINATION_ID),]

sub2=data.frame(COMBINATION_ID=m$COMBINATION_ID, CELL_LINE=m$CELL_LINE, prediction=pnorm( o$par$ytest/sqrt(o$par$sigma_sq ) ) )

sub2mat=dcast(sub2, COMBINATION_ID ~ CELL_LINE, value.var="prediction")

rownames(sub2mat)=sub2mat$COMBINATION_ID
sub2mat$COMBINATION_ID=NULL

dir.create("sub2")
setwd("sub2")
write.csv(sub2mat, "confidence_matrix.csv", quote=F)

sub2bin =sub2mat>.9
class(sub2bin)="numeric"
write.csv(sub2bin, "synergy_matrix.csv", quote=F)

sub1=data.frame(CELL_LINE=test$CELL_LINE, COMBINATION_ID=test$COMBINATION_ID, PREDICTION=o$par$ytest)
dir.create("sub1")
write.csv(sub1,file="sub1/prediction.csv",row.names = F, quote=F)

system("zip sub2.zip *.csv")
setwd("..")
sub1conf=data.frame(COMBINATION_ID=a$COMBINATION_ID, CONFIDENCE=1-pnorm( -abs(a$mean) / sqrt(a$err) )*2 )
write.csv(sub1conf,file="sub1/combination_priority.csv",row.names = F, quote=F)
