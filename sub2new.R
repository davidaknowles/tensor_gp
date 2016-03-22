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
    setup="sub2final"
    sub_challenge="2"
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
  
  dat$Ntest=nrow(m)
  dat$cellLinesTest=as.numeric(m$CELL_LINE)
  dat$drugATest=as.numeric(m$COMPOUND_A)
  dat$drugBTest=as.numeric(m$COMPOUND_B)

stop()

cat("Compiling...\n")
  sm=stan_model("comb_therapy_models/gp_multitask_mkl.stan")

cat("Running...\n")
  load("cached_results/sub2_sub2final_tissue1_seed1_iter30.RData")

o_mat=optimizing(sm, data=dat, verbose=T, init=o$par[1:8], as_vector=F, iter=1)

cat("Saving...\n")
  save(o_mat, file="cached_results/sub2_sub2final_all.RData")

cat("Done\n")
