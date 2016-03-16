require(rstan)
require(doMC)
registerDoMC(3)
source("load_data.R")
source("load_sub1final_data.R")

sub1b=F

sm=stan_model("comb_therapy_models/gp_multitask_mkl.stan")

cls=levels(train$CELL_LINE)

if (sub1b) {
  dist=dist[c("cnv","mut")]
}
# only allowed CNV and MUT
sqDist=lapply(,function(g) g[cls,cls]^2)

drugs=levels(train$COMPOUND_A)

pathways=read.csv("DREAM CHALLENGE TABLE DRUGS-TARGETS LOUKIA.csv",row.names = 1,check.names = F)
pathways=pathways[,colSums(pathways)>1]
ldrugs=do.call(rbind,strsplit(rownames(pathways),"-"))[,2]
sort(drugs[ !drugs %in% ldrugs])
sort( ldrugs[ ! ldrugs %in% drugs ] )

rownames(pathways)=ldrugs
pathways=pathways[drugs,]
pathways=as.matrix(pathways)
#a=pathways %*% t(pathways)

require(sna)
graph_dist=geodist(pathways,count.paths = F, inf.replace = 200)$gdist[1:nrow(pathways),1:nrow(pathways)]/2
graph_dist=graph_dist/median(graph_dist[upper.tri(graph_dist)])

# not allowed monotherapy data
sqDist_dr=list( graph_dist^2  ) 

dat=list(N=nrow(train), Ntest=nrow(test), y=train$SYNERGY_SCORE, C=length(cls), D=length(levels(train$COMPOUND_A)), P=length(sqDist), sqDist_cl=sqDist, P_dr=length(sqDist_dr), sqDist_dr=sqDist_dr, cellLines=as.numeric(train$CELL_LINE), cellLinesTest=as.numeric(test$CELL_LINE), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), drugATest=as.numeric(test$COMPOUND_A), drugBTest=as.numeric(test$COMPOUND_B) )


reruns = foreach(i=1:6) %dopar% { optimizing(sm, data=dat,verbose=T,as_vector=F) }

save(reruns, file="sub1_partB_final.RData")

likelihoods=foreach(r=reruns, .combine = c) %do% r$value
barplot(-likelihoods)

o=reruns[[ which.max(likelihoods) ]]

# o$value
require(gridExtra)
do.call(grid.arrange , c(foreach(r=reruns[order(likelihoods)]) %do% {
  ggplot(data.frame(x=names(sqDist), y=sqrt(r$par$eta_sq_cl)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("Likelihood: ",format(r$value, digits = 3))) + theme_bw(base_size=14) + ylab("importance") + xlab("") + ylim(0,2.3) } , nrow=1))

#names(o$par$eta_sq_dr)=c("pathways","mono_therapy")
#barplot(sqrt(o$par$eta_sq_dr),ylab="importance")

N=nrow(train)
o$par$sigma_sq
SN=o$par$Sigma_no_noise

ytrain = SN %*% ( solve(SN + diag(N) * o$par$sigma_sq, dat$y - o$par$mu) ) + o$par$mu
qplot(dat$y, ytrain)
train$err= (ytrain-dat$y)^2

require(dplyr)
errs=train %>% group_by(COMBINATION_ID) %>% summarize(err=mean(err), mean=mean(SYNERGY_SCORE))
setDF(errs)
a=as.data.frame(errs[1:167,])

plot( a$mean, sqrt(a$err ), xlab="mean synergy", ylab="training RMSE", pch=16 ) 
#text( a$mean, sqrt(a$err ), a$COMBINATION_ID)
a$cv=sqrt(a$err)/ abs(a$mean) 
a=a[order(a$cv),]

b=a[c(1:20,nrow(a)-20+(1:20)),]
b$COMBINATION_ID=factor(b$COMBINATION_ID,b$COMBINATION_ID)
ggplot(b, aes(x=COMBINATION_ID, y=cv) ) + geom_bar(stat = "identity")+coord_flip() + scale_y_log10()

str(o)
sub1=data.frame(CELL_LINE=test$CELL_LINE, COMBINATION_ID=test$COMBINATION_ID, PREDICTION=o$par$ytest)

resdir="sub1_partB_final"
setwd(resdir)
write.csv(sub1,file="prediction.csv",row.names = F, quote=F)

sub1conf=data.frame(COMBINATION_ID=a$COMBINATION_ID, CONFIDENCE=1-pnorm( -abs(a$mean) / sqrt(a$err) )*2 )
write.csv(sub1conf,file="combination_priority.csv",row.names = F, quote=F)

system("zip sub1b.zip *.csv")
setwd("..")
