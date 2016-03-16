source("mono_therapy_models/mono_model.R")
source("utils.R")

require(doMC)
#registerDoMC(if (detectCores()==8) 7 else detectCores())
registerDoMC(2)

cls=levels(mono$CELL_LINE)
dat=list(K=10,D=length(levels(mono$COMPOUND)), C=length(levels(mono$CELL_LINE)), N=nrow(mono), drug=as.numeric(mono$COMPOUND), cellLine=as.numeric(mono$CELL_LINE), y=mono$Einf, Ntest=nrow(monotest), drugtest=as.numeric(monotest$COMPOUND), cellLinetest=as.numeric(monotest$CELL_LINE), sqDist=gedist[cls,cls]^2)

#-------------------
#mf =stan_model("mono_therapy_models/matrix_factorization_gp_alt.stan")

mf =stan_model("mono_therapy_models/matrix_factorization_gp_gp_tau_local.stan")
res=foreach (iterat=1:10) %dopar% {
  v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=1200, output_samples=30, seed=iterat)
  getMeans(v)
}
#save(file="gpgp_K10_restarts.RData",res)
cat("Done\n")


load("gpgp_K10_restarts.RData")
rmsesReps=unlist(lapply(res, function(g) rmseOnMonoTest(g$testMean)))
g=res[[which.min(rmsesReps)]]
drugV=g$drugV
rownames(drugV)=levels(mono$COMPOUND)
apply(drugV,2,sd)
require(irlba)
s=irlba(drugV)
ggplot(data.frame(pc1=s$u[,1],pc2=s$u[,2],name=levels(mono$COMPOUND)),aes(pc1,pc2,label=name))+geom_text()

boxplot(rmses,ylab="Test set RMSE")

iterations=seq(100,1200,by=100)
rmses=unlist(foreach (iterat=iterations) %dopar% {
  v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=iterat, output_samples=30, seed=1)
  rmseOnMonoTest(getMeans(v)$testMean)
})
save(file="gpgp_K10_progression.RData",iterations,rmses)

dat$K=10
v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=1000, output_samples=30, seed=1)
g=getMeans(v)
rmseOnMonoTest(getMeans(v)$test_synergy)
save(g,file="gpgp1000.RData")
cat("DONE\n")

#---------------- NO CELL LINE DATA -----------------------
mf =stan_model("mono_therapy_models/matrix_factorization.stan")
numFact=floor(seq(1,10.9,by=.1))
rmses=unlist(foreach (k=numFact) %dopar% {
  dat$K=k
  v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=800, output_samples=100)
  rmseOnMonoTest(getMeans(v)$test_synergy)
})
names(rmses)=numFact
save(rmses,file="mono_therapy_models/mf_no_cl_data.RData")

iterations=seq(100,700,by=100)
rmses=unlist(foreach (iterat=iterations) %dopar% {
  v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=1000, output_samples=30, seed=1)
  rmseOnMonoTest(getMeans(v)$testMean)
})
qplot( c( iterations, rep(1200,10)) ,c( rmses, rmsesReps))+xlab("iterations")+ylab("Test set RMSE")+theme_bw(base_size = 18)+geom_line(size=1)+geom_point(size=3)


numFact=floor(seq(1,10.9,by=.1))
rmses=unlist(foreach (k=numFact) %dopar% {
  dat$K=k
  v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=800, output_samples=100)
  rmseOnMonoTest(getMeans(v)$test_synergy)
})
names(rmses)=numFact
qplot(as.factor(numFact), rmses)+geom_boxplot()+theme_bw(base_size = 18)+xlab("# factors")+ylab("Test set RMSE")

eightFact=foreach (k=1:10) %dopar% {
  dat$K=8
  v=rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=800, output_samples=100)
  getMeans(v)
}
eightFactRmses=sapply(eightFact,function(g) rmseOnMonoTest(g$test_synergy))
best=eightFact[[which.min(eightFactRmses)]]
drugF=best$drugV
rownames(drugF)=levels(mono$COMPOUND)
plot(hclust(dist(drugF)))
s=svd(drugF)
plot(s$d^2/sum(s$d^2))

plot(s$u[,1],s$u[,2],col="white")
text(s$u[,1],s$u[,2],levels(mono$COMPOUND))
apply(best$cellV,1,sd)

# very slow!
mf =stan_model("mono_therapy_models/matrix_factorization_gp_gp.stan")

mf =stan_model("mono_therapy_models/matrix_factorization_gp_nn.stan")
v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=300, output_samples=100)
rmseOnMonoTest(getMeans(v)$testMean)

numFact=floor(seq(1,10.9,by=.1))
dat$K2=5
y=scale(dat$y,scale=F)
dat$y=as.numeric(y)
rmses=unlist(foreach (k=numFact) %dopar% {
  dat$K=k
  v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=800, output_samples=100)
  rmseOnMonoTest(getMeans(v)$testMean)
})
names(rmses)=numFact
qplot(as.factor(numFact), rmses)+geom_boxplot()+theme_bw(base_size = 18)+xlab("# factors")+ylab("Test set RMSE")

mf=stan_model("mono_therapy_models/matrix_factorization_l2.stan")
numFact=floor(seq(1,10.9,by=.1))
dat$sqDist=NULL
ges=scale(t(ge))
xx=ges %*% t(ges)
xx=xx/mean(diag(xx))
dat$xx=xx[levels(mono$CELL_LINE),levels(mono$CELL_LINE)]
rmses=unlist(foreach (k=numFact) %dopar% {
  dat$K=k
  tryCatch( {
    v = rstan:::vb(mf, data=dat, eta_adagrad=0.01, iter=800, output_samples=100)
    rmseOnMonoTest(getMeans(v)$test_synergy)
  }, error=function(g) NA )
})
names(rmses)=numFact
qplot(as.factor(numFact), rmses)+geom_boxplot()+theme_bw(base_size = 18)+xlab("# factors")+ylab("Test set RMSE")

#qplot(as.factor(numFact), rmses)+geom_boxplot()+theme_bw(base_size = 18)+xlab("# factors")+ylab("Test set RMSE")

load("")
