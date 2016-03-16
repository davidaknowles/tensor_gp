require(dplyr)
require(Matrix)

source("convex_pca.R")

convert_to_matrices=function(input) {
  means=input %>% group_by(CELL_LINE,COMPOUND) %>% summarize(mean=mean(Einf), n=length(Einf)) %>% as.data.frame
  list( means=as.matrix(sparseMatrix(i=as.numeric(means$CELL_LINE),j=as.numeric(means$COMPOUND),x=means$mean)),
weights=as.matrix(sparseMatrix(i=as.numeric(means$CELL_LINE),j=as.numeric(means$COMPOUND),x=means$n)) )
}

training=convert_to_matrices(mono)
testing=convert_to_matrices(monotest)

best_r=cv.convexPCA( training$means, weights=training$weights )
co=convexPCA( training$means, best_r, weights=training$weights, Ytest=testing$means, test_weights = testing$weights, verbose=T, rmseTol=1e-8)

comb_data=rbind(mono,monotest)
combined=convert_to_matrices(comb_data)
#combined$means=combined$means-mean(comb_data$Einf)
best_r=cv.convexPCA( combined$means, weights=combined$weights )
co=convexPCA( combined$means, best_r, weights=combined$weights,  verbose=T, rmseTol=1e-8)
#co=co+mean(comb_data$Einf)
dimnames(co)=list( levels(mono$CELL_LINE), levels(mono$COMPOUND) )
require(gplots)
heatmap.2(t(co),trace="none",symbreaks=F, key.title="", margins=c(5,8))

top_two_pcs=irlba(co,2,2)
plot(top_two_pcs$v[,1], top_two_pcs$v[,2], col="white", xlab="PC 1", ylab="PC 2")
text(top_two_pcs$v[,1], top_two_pcs$v[,2], colnames(co))

#write.table(co,"pca_imputed_mono.txt",quote=F)

mono_dist=as.matrix(dist(t(co)))
mono_dist=mono_dist/median(mono_dist[upper.tri(mono_dist)])
