
views=c("cnv", "gex", "mut", "methyl", "tissue")
dist=foreach(view=views) %do% {
  a=read.csv(paste0("processed_data/",view,"_dist.csv"),row.names=1,check.names=F)
  a=a/median(a[upper.tri(a)])
  as.matrix(a)
}
names(dist)=views
dist=lapply(dist, function(g) g[rownames(dist$gex),colnames(dist$gex)])