
require(irlba)
pc=irlba(as.matrix(ge), 30, 30)$v
rownames(pc)=colnames(ge)

nn=stan_model("nn_real.stan")

as.matrix(model.matrix(~COMPOUND,data=mono)) -> x
x=x[,2:ncol(x)]
X=cbind(x,pc[mono$CELL_LINE,])

as.matrix(model.matrix(~COMPOUND,data=monotest)) -> xt
xt=xt[,2:ncol(xt)]
Xt=cbind(xt,pc[monotest$CELL_LINE,])

y=scale(mono$Einf)
atty=attributes(y)
dat=list(num_nodes=10, num_middle_layers=4, d=ncol(X), N=nrow(X), Nt=nrow(Xt), X=X, y=y, Xt=Xt)
#s <- sampling(nn, data = dat, iter = 1000, chains = 4, cores=4)
v <- rstan:::vb(nn, data = dat, iter=1000,  output_samples=100)
preds=getMeans(v)$predictions
atty=attributes(y)
preds=preds*atty[["scaled:scale"]]+atty[["scaled:center"]]
rmseOnMonoTest(preds) # 37.2

save(s,file="nn_mono.RData")

load("")

stop()

load("nn_mono.RData")
samp=extract(s)
preds=colMeans(samp$predictions)
rmseOnMonoTest(preds) # 57.7