require(rstan)
require(doMC)

if (!interactive()) {
    args=commandArgs(trailingOnly = T)
    run=as.logical(as.numeric(args[1])) # 0 or 1, 1=run the model, 0=just load cached results
    setup=args[2] # lb or final, or sub2
    sub_challenge=args[3] # A, B or 2
    registerDoMC( as.numeric(args[4] )) # number of cores
    use_tissue=as.logical(as.numeric(args[5])) # 0/1, use tissue similarity? 
} else {
    run=F
    setup="sub2final"
    sub_challenge="2"
    use_tissue=T
}
source("load_cell_line_data.R")
source("load_response_data.R")

dat=switch(setup, 
  lb=load_data( "ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv"), 
  lb2=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch2_LB.csv"),"ch1_LB.csv"), 
  final=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv"),"ch1_leaderBoard_monoTherapy.csv"),
  final2=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv","ch2_LB.csv"),"ch1_leaderBoard_monoTherapy.csv"),
  sub2=load_data( "ch1_train_combination_and_monoTherapy.csv","ch2_LB.csv"),
  sub2final=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv","ch2_LB.csv"),"ch2_test_monoTherapy.csv"),
)

train=dat$train
test=dat$test

cls=levels(train$CELL_LINE)

if (!use_tissue)
  dist=dist[ names(dist) != "tissue" ]

if (sub_challenge=="B") 
  dist=dist[ ! (names(dist) %in% c("gex","methyl") ) ]

sqDist=lapply(dist,function(g) g[cls,cls]^2)

drugs=levels(train$COMPOUND_A)

pathways=read.csv("processed_data/drug_targets.csv",row.names = 1,check.names = F)
pathways=pathways[,colSums(pathways)>1]
ldrugs=do.call(rbind,strsplit(rownames(pathways),"-"))[,2]
sort(drugs[ !drugs %in% ldrugs])
sort( ldrugs[ ! ldrugs %in% drugs ] )

rownames(pathways)=ldrugs
pathways=pathways[drugs,]
pathways=as.matrix(pathways)

require(sna)
graph_dist=geodist(pathways,count.paths = F, inf.replace = 200)$gdist[1:nrow(pathways),1:nrow(pathways)]/2
graph_dist=graph_dist/median(graph_dist[upper.tri(graph_dist)])

sqDist_dr=list( pathway=graph_dist^2 ) # sdm

if (sub_challenge=="A" || sub_challenge=="2") {
  co=read.table("processed_data/pca_imputed_mono.txt",header=T)
  mono_dist=as.matrix(dist(t(co)))
  mono_dist=mono_dist/median(mono_dist[upper.tri(mono_dist)])
  sqDist_dr$mono=mono_dist[drugs,drugs]^2
}

dat=list(N=nrow(train), Ntest=nrow(test), y=train$SYNERGY_SCORE, C=length(cls), D=length(levels(train$COMPOUND_A)), P=length(sqDist), sqDist_cl=sqDist, P_dr=length(sqDist_dr), sqDist_dr=sqDist_dr, cellLines=as.numeric(train$CELL_LINE), cellLinesTest=as.numeric(test$CELL_LINE), drugA=as.numeric(train$COMPOUND_A), drugB=as.numeric(train$COMPOUND_B), drugATest=as.numeric(test$COMPOUND_A), drugBTest=as.numeric(test$COMPOUND_B) )

if (run) {

    sm=stan_model("comb_therapy_models/gp_multitask_mkl.stan")

  foreach(i=1:10) %dopar% {
      
    resfile=paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",i,".RData")
    if (file.exists(resfile)) return(NULL)
    set.seed(i)
    o=optimizing(sm, data=dat,verbose=T,as_vector=F) 
    save(o, file=resfile)
      
  }
  cat("Done!")
  stop()
} else {
  reruns = foreach(i=1:10) %dopar% { 
    load(paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",i,".RData"))
    o
  }
}

likelihoods=foreach(r=reruns, .combine = c) %do% r$value
barplot(-likelihoods)

source("utils.R")
scores_df=foreach(i=seq_along(reruns), .combine = rbind) %do% { reruns[[i]]$score=get_score(reruns[[i]]$par$ytest, test)
                                                         c( reruns[[i]]$score, attr(reruns[[i]]$score, "se")) }
scores=scores_df[,1]
scores[is.na(scores)]=0

o=reruns[[ which.max(likelihoods) ]]
y=c(with_tissue, o_score)
se=c(attr(with_tissue,"se"), attr(o_score,"se"))
qplot( c("with tissue","without"), y, geom="blank" ) + geom_bar(stat="identity", position="dodge") + theme_bw(base_size = 16) + xlab("") + ylab("Score") + geom_errorbar( aes(ymax=y+se, ymin=y-se), width=.6 )

# o$value
require(gridExtra)
do.call(grid.arrange , c(foreach(r=reruns[order(likelihoods)]) %do% { ggplot(data.frame(x=names(sqDist)[1:4], y=sqrt(r$par$eta_sq_cl)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("L: ",format(r$value, digits = 3)," S:",r$score)) + theme_bw(base_size=14) + ylab("importance") + xlab("") + ylim(0,2.3) } , nrow=2))

ggplot(data.frame(x=names(sqDist), y=sqrt(o$par$eta_sq_cl)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("Likelihood: ",format(r$value, digits = 3))) + theme_bw(base_size=16) + ylab("importance") + xlab("") + ylim(0,2.3)

do.call(grid.arrange , c(foreach(r=reruns[order(likelihoods)]) %do% { ggplot(data.frame(x=c("pathways","mono_therapy"), y=sqrt(r$par$eta_sq_dr)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("L: ",format(r$value, digits = 3)," S:",r$score))+ theme_bw(base_size=14) + ylab("importance") + xlab("") + ylim(0,2.3) } , nrow=2))

ggplot(data.frame(x=c("pathways","mono_therapy"), y=sqrt(o$par$eta_sq_dr)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("Likelihood: ",format(r$value, digits = 3))) + theme_bw(base_size=16) + ylab("importance") + xlab("") + ylim(0,2.3)

N=nrow(train)
o$par$sigma_sq
SN=o$par$Sigma_no_noise

ytrain = SN %*% ( solve(SN + diag(N) * o$par$sigma_sq, dat$y - o$par$mu) ) + o$par$mu
qplot(dat$y, ytrain) + theme_bw(base_size = 16) + xlab("True synergy") + ylab("Predicted") + geom_abline(intercept = 0, slope=1)
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

ggplot(data.frame(likelihoods, scores), aes(likelihoods, scores)) +geom_point(size=3) + theme_bw(base_size = 16) + ylab("Score") + xlab("Log likelihood")

o_score=get_score(o$par$ytest, test)
df=attr(o_score, "df")
dfsub=df[df$n>2,]
ggplot(dfsub, aes(as.factor(n), pearson)) + geom_boxplot(outlier.shape = NA) + xlab("# cell lines in combination") + ylab("Pearson correlation") + geom_point(position = position_jitter(width = .5), size=3, alpha=.5) + theme_bw(base_size = 15)

str(o)
sub1=data.frame(CELL_LINE=test$CELL_LINE, COMBINATION_ID=test$COMBINATION_ID, PREDICTION=o$par$ytest)

resdir=paste0("sub1_part",sub_challenge,"_final")
dir.create(resdir)
setwd(resdir)
write.csv(sub1,file="prediction.csv",row.names = F, quote=F)

sub1conf=data.frame(COMBINATION_ID=a$COMBINATION_ID, CONFIDENCE=1-pnorm( -abs(a$mean) / sqrt(a$err) )*2 )
write.csv(sub1conf,file="combination_priority.csv",row.names = F, quote=F)

system(paste0("zip sub1",sub_challenge,".zip *.csv"))
setwd("..")




if (subchallenge=="sub2") {
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
}

