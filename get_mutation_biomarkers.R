require(rstan)
require(doMC)
#sm=stan_model("comb_therapy_models/gp_multitask_mkl_importance.stan")

setup="final2"
sub_challenge="B"
use_tissue=T
iterations=30


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

reruns = foreach(i=1:10) %dopar% { 
  fn=paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",i,"_iter",iterations,".RData")
  if (!file.exists(fn)) return(NULL)
  load(fn)
  o
}
likelihoods=foreach(r=reruns, .combine = c) %do% r$value
plot(likelihoods)
o=reruns[[ which.max(likelihoods) ]]
rm("reruns","r"); gc()
#fn=paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",1,"_iter",iterations,".RData")
#load(fn)

SN=o$par$Sigma_no_noise
N=nrow(SN)

names(o$par$eta_sq_cl)=names(sqDist)
names(o$par$inv_rho_sq_cl)=names(sqDist)

view_to_investigate="mut"

sqDist=sqDist[ ! (names(sqDist) %in% view_to_investigate ) ]

dat=list(N=nrow(train), y=train$SYNERGY_SCORE, C=length(cls), D=length(levels(train$COMPOUND_A)), P=length(sqDist), sqDist_cl=sqDist, P_dr=length(sqDist_dr), sqDist_dr=sqDist_dr, cellLines=as.integer(train$CELL_LINE),  drugA=as.integer(train$COMPOUND_A), drugB=as.integer(train$COMPOUND_B) )

dat$eta_sq_cl=o$par$eta_sq_cl[names(sqDist)]
dat$inv_rho_sq_cl=o$par$inv_rho_sq_cl[names(sqDist)]
dat=c(dat, o$par[c("sigma_sq_cl","eta_sq_dr","inv_rho_sq_dr","sigma_sq_dr")])

# TODO: make this condition on mut vs. ge
dat$eta_sq_cl_mut=o$par$eta_sq_cl[ view_to_investigate ]
dat$inv_rho_sq_mut=o$par$inv_rho_sq_cl[ view_to_investigate ]

dat$unmixed_y = solve(SN + diag(N) * o$par$sigma_sq, dat$y - o$par$mu) 

require(data.table)

a=fread("~/Box Sync/astrazeneca_dream/code/Sanger_molecular_data/mutations.csv")
setDF(a)
a$cell_line_name=factor(a$cell_line_name)
a$Gene.name=factor(a$Gene.name)
require(Matrix)
sm=sparseMatrix(i=as.numeric(a$cell_line_name), j=as.numeric(a$Gene.name), x=rep(1,nrow(a)), dimnames = list(levels(a$cell_line_name),levels(a$Gene.name)))
mut=t(as.matrix(sm))>0
rm(a); gc()

require(prabclus)
d=jaccard(mut)
dimnames(d)=list(rownames(sm),rownames(sm))

median_d=median(d[upper.tri(d)])

dat$inv_rho_sq_mut = dat$inv_rho_sq_mut / median_d^2

mode(mut)="numeric"
dat$mut=t(mut[,cls])
dat$num_genes=ncol(dat$mut)

combs_to_test=read.table("biomarker_to_predict.txt", header=F, stringsAsFactors = F)$V1
comb=combs_to_test[1]

csv_lines=foreach(comb=combs_to_test) %do% {
  wh=which(train$COMBINATION_ID==comb)
  dat$drugATest=as.integer(train[wh[1],"COMPOUND_A"])
  dat$drugBTest=as.integer(train[wh[1],"COMPOUND_B"])
  
  #registerDoMC(7)
  g=foreach(i=wh, .combine = cbind) %do% {
    # testing one combination
    dat$cellLinesTest=as.integer(train[i,"CELL_LINE"])
    s=stan("comb_therapy_models/gp_multitask_mkl_mut_importance.stan", data=dat, chains=0)
    grad_log_prob(s, mut[,as.integer(train[i,"CELL_LINE"])])
  }

  pv=2*foreach(gene_index=seq_len(nrow(g)), .combine = c) %dopar% 
    pnorm( abs(mean( g[gene_index,] ) / sd(g[gene_index,])), lower.tail = F )
  
  gi=order(pv)[3]
  
  cl=as.character(train[wh, "CELL_LINE"])
  qplot( mut[gi,cl], train[wh, "SYNERGY_SCORE"], label=cl, geom="text" ) + xlab(paste(rownames(mut)[gi])) + ylab(paste(comb,"SYNGERY")) + theme_bw(base_size=20)
  
  top_hits=order(pv)[1:10]
  
  direction=foreach(gene_index=top_hits, .combine = c) %do% mean( g[gene_index,] )
  
  df=data.frame( gene=rownames(mut)[ top_hits ], type="mutation", feat=1, direct=ifelse( direction > 0, 1, -1))
  paste(as.character( t(as.matrix(df)) ), collapse = ",")
} 

rownames(mut)[ order(pv)[1:10] ]