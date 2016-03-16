
# ---------- GPGP
dat=list(K=10, D=length(levels(mono$COMPOUND)), C=length(levels(mono$CELL_LINE)), N=nrow(mono), drug=as.numeric(mono$COMPOUND), cellLine=as.numeric(mono$CELL_LINE), y=mono$Einf, Ntest=nrow(monotest), drugtest=as.numeric(monotest$COMPOUND), cellLinetest=as.numeric(monotest$CELL_LINE), sqDist=gedist[cls,cls]^2)

mf=stan_model("matrix_factorization_gp_gp.stan")
#s=sampling(mf, data=dat, chains=7, cores=7, iter=1000, verbose=T)
o=optimizing(mf, data=dat, verbose=T, iter=8, as_vector=F)
s=sampling(mf, data=dat, chains=1, cores=1, iter=100, verbose=T)
samps=extract(s)

preds=colMeans(samps$test_synergy)
rmseOnMonoTest(preds)

