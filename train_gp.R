require(rstan)
require(doMC)

if (!interactive()) {
    args=commandArgs(trailingOnly = T)
    if (length(args)<6)
        stop("Usage is: Rscript train_gp.R <run> <setup> <sub> <cores> <usetissue> <max_its>\n where \n run: whether to train the model(s) or analyze the results. \n setup: lb, lb2, final, final2, sub2 or sub2final (see code) \n sub: A, B or 2. \n cores: number of cores to use. \n usetissue: 0/1")
    run=as.logical(as.numeric(args[1])) # 0 or 1, 1=run the model, 0=just load cached results
    setup=args[2] # lb or final, or sub2
    sub_challenge=args[3] # A, B or 2
    cores=as.numeric(args[4] )
    if (cores>1) registerDoMC(cores) # number of cores
    use_tissue=as.logical(as.numeric(args[5])) # 0/1, use tissue similarity? 
    iterations=as.numeric(args[6])
} else {
    run=T
    setup="final2"
    sub_challenge="A"
    use_tissue=T
    iterations=30
}

source("load_cell_line_data.R")
source("load_response_data.R")

dat=switch(setup, 
  lb=load_data( "ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv"), # original leaderboard
  lb2=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch2_LB.csv"),"ch1_LB.csv"), # leaderboard also using training from Challenge 2
  final=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv"),"ch1_test_monoTherapy.csv"), #  final
  final2=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv","ch2_LB.csv"),"ch1_test_monoTherapy.csv"), # final using Ch 2 data in addition
  sub2=load_data( "ch1_train_combination_and_monoTherapy.csv","ch2_LB.csv"), # Challenge 2 leaderboard
  sub2final=load_data( c("ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv","ch2_LB.csv"),"ch2_test_monoTherapy.csv"), # Challenge 2 final
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
      
    resfile=paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",i,"_iter",iterations,".RData")
    if (file.exists(resfile)) return(NULL)
    set.seed(i)
    st=system.time({ o=optimizing(sm, data=dat, verbose=T, init=init, as_vector=F, iter=iterations) })[1]
    attr(o,"time")=st
    save(o, file=resfile)
  }
  cat("Done!")
  stop()
} else {
  reruns = foreach(i=1:10) %dopar% { 
    fn=paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",i,"_iter",iterations,".RData")
    if (!file.exists(fn)) return(NULL)
    load(fn)
    o
  }
}

likelihoods=foreach(r=reruns, .combine = c) %do% r$value
plot(likelihoods)
o=reruns[[ which.max(likelihoods) ]]

if (0) { # scoring, only for leaderboard, not final run of course
  source("utils.R")
  scores_df=foreach(i=seq_along(reruns), .combine = rbind) %do% { reruns[[i]]$score=get_score(reruns[[i]]$par$ytest, test)
                                                           c( reruns[[i]]$score, attr(reruns[[i]]$score, "se")) }
  scores=scores_df[,1]
  scores[is.na(scores)]=0
}

if (0) { # plot importance of different distance matrices for each run
  require(gridExtra)
  do.call(grid.arrange , c(foreach(r=reruns[order(likelihoods)]) %do% { ggplot(data.frame(x=names(sqDist), y=sqrt(r$par$eta_sq_cl)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("L: ",format(r$value, digits = 3)," S:",r$score)) + theme_bw(base_size=14) + ylab("importance") + xlab("") + ylim(0,2.3) } , nrow=2))
  
  ggplot(data.frame(x=names(sqDist), y=sqrt(o$par$eta_sq_cl)), aes(x,y)) + geom_bar(stat="identity") + theme_bw(base_size=16) + ylab("importance") + xlab("") 
  
  do.call(grid.arrange , c(foreach(r=reruns[order(likelihoods)]) %do% { ggplot(data.frame(x=c("pathways","mono_therapy"), y=sqrt(r$par$eta_sq_dr)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("L: ",format(r$value, digits = 3)," S:",r$score))+ theme_bw(base_size=14) + ylab("importance") + xlab("") + ylim(0,2.3) } , nrow=2))
  
  ggplot(data.frame(x=c("pathways","mono_therapy"), y=sqrt(o$par$eta_sq_dr)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("Likelihood: ",format(r$value, digits = 3))) + theme_bw(base_size=16) + ylab("importance") + xlab("") 
}

N=nrow(train)
SN=o$par$Sigma_no_noise

ytrain = SN %*% ( solve(SN + diag(N) * o$par$sigma_sq, dat$y - o$par$mu) ) + o$par$mu
train$SE= (ytrain-dat$y)^2

require(dplyr)
training_errors=train %>% group_by(COMBINATION_ID) %>% summarize(n=length(SE), RMSE=sqrt(sum(SE)/(n-1)),  mean=mean(SYNERGY_SCORE), var=var(SYNERGY_SCORE))
setDF(training_errors)

if (0) { # just some plotting

  qplot(dat$y, ytrain) + theme_bw(base_size = 16) + xlab("True synergy") + ylab("Predicted") + geom_abline(intercept = 0, slope=1)
  
  plot( training_errors$mean, training_errors$RMSE, xlab="mean synergy", ylab="training RMSE", pch=16 ) 
  
  training_errors$cv=training_errors$RMSE/ abs(training_errors$mean) 
  training_errors=training_errors[order(training_errors$cv),]
  
  b=training_errors[c(1:20,nrow(training_errors)-20+(1:20)),]
  b$COMBINATION_ID=factor(b$COMBINATION_ID,b$COMBINATION_ID)
  ggplot(b, aes(x=COMBINATION_ID, y=cv) ) + geom_bar(stat = "identity")+coord_flip() + scale_y_log10()
  
  ggplot(data.frame(likelihoods, scores), aes(likelihoods, scores)) +geom_point(size=3) + theme_bw(base_size = 16) + ylab("Score") + xlab("Log likelihood")
  
  o_score=get_score(o$par$ytest, test)
  df=attr(o_score, "df")
  dfsub=df[df$n>2,]
  ggplot(dfsub, aes(as.factor(n), pearson)) + geom_boxplot(outlier.shape = NA) + xlab("# cell lines in combination") + ylab("Pearson correlation") + geom_point(position = position_jitter(width = .5), size=3, alpha=.5) + theme_bw(base_size = 15)
}

if (sub_challenge %in% c("A","B")) {

  predictions=data.frame(CELL_LINE=test$CELL_LINE, COMBINATION_ID=test$COMBINATION_ID, PREDICTION=o$par$ytest)
  
  resdir=paste0("sub1_part",sub_challenge,"_",setup)
  dir.create(resdir)
  setwd(resdir)
  write.csv(predictions,file="prediction.csv",row.names = F, quote=F)
  
  training_errors=training_errors[ as.character(training_errors$COMBINATION_ID) %in% as.character(test$COMBINATION_ID), ]
  
  #sub1conf=data.frame(COMBINATION_ID=training_errors$COMBINATION_ID, CONFIDENCE=1-pnorm( -abs(training_errors$mean) / training_errors$RMSE )*2 )
  sub1conf=data.frame(COMBINATION_ID=training_errors$COMBINATION_ID, CONFIDENCE=rank( training_errors$var / (training_errors$RMSE^2) )/nrow(training_errors) )
  write.csv(sub1conf,file="combination_priority.csv",row.names = F, quote=F)
  
  system(paste0("zip sub1",sub_challenge,".zip *.csv"))
  setwd("..")
}
