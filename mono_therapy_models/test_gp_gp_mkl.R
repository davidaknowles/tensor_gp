source("mono_therapy_models/mono_model.R")
source("utils.R")

require(doMC)
#registerDoMC(if (detectCores()==8) 7 else detectCores())
registerDoMC(2)

cls=levels(mono$CELL_LINE)
dat=list(P=length(dist),K=10,D=length(levels(mono$COMPOUND)), C=length(levels(mono$CELL_LINE)), N=nrow(mono), drug=as.numeric(mono$COMPOUND), cellLine=as.numeric(mono$CELL_LINE), y=(mono$Einf-50)/10, Ntest=nrow(monotest), drugtest=as.numeric(monotest$COMPOUND), cellLinetest=as.numeric(monotest$CELL_LINE), sqDist=lapply(dist,function(g) g[cls,cls]^2))

#-------------------
#mf =stan_model("mono_therapy_models/matrix_factorization_gp_alt.stan")

mf =stan_model("mono_therapy_models/matrix_factorization_gp_gp_mkl.stan")

v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=1000, output_samples=30, seed=1)
tm=getMeans(v)
rmseOnMonoTest(tm$testMean*10+50)
save(file="gpgp_mkl_2.RData",tm)
cat("Done\n")
