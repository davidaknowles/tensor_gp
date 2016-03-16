
source("load_data.R")

bin=do.call(rbind,foreach(i=1:nrow(ge)) %dopar% {
  g=as.numeric(ge[i,])
  km=kmeans( g, centers=c(min(g),max(g)))
  if (km$centers[2] > km$centers[1]) km$cluster-1 else 1-km$cluster
})

cors=unlist(foreach(i=1:nrow(ge)) %do% cor( as.numeric(bin[i,]),as.numeric(ge[i,])))

d=as.matrix(dist(t(bin)))
d=d/median(d[upper.tri(d)])
sd(d[upper.tri(d)])
hist(d[upper.tri(d)])

j=1-jaccard(t(bin))
j=j/median(j[upper.tri(j)])
sd(j[upper.tri(j)])
hist(j[upper.tri(j)])
