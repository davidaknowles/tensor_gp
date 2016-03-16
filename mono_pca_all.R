
require(dplyr)
require(Matrix)
source("convex_pca.R")

drug_dir="Drug_synergy_data/"
all_mono=foreach (f=list.files(drug_dir,"_monoTherapy.csv"), .combine = rbind) %do% getMono(read.csv(paste0(drug_dir,f),stringsAsFactors = F))

all_mono=all_mono[!duplicated(all_mono),]

all_mono$CELL_LINE=as.factor(all_mono$CELL_LINE)
all_mono$COMPOUND=as.factor(all_mono$COMPOUND)

convert_to_matrices=function(input) {
  means=input %>% group_by(CELL_LINE,COMPOUND) %>% summarize(mean=mean(Einf), n=length(Einf)) %>% as.data.frame
  list( means=as.matrix(sparseMatrix(i=as.numeric(means$CELL_LINE),j=as.numeric(means$COMPOUND),x=means$mean)),
weights=as.matrix(sparseMatrix(i=as.numeric(means$CELL_LINE),j=as.numeric(means$COMPOUND),x=means$n)) )
}

training=convert_to_matrices(all_mono)

best_r=cv.convexPCA( training$means, weights=training$weights )
co=convexPCA( training$means, best_r, weights=training$weights, verbose=T, rmseTol=1e-8, its = 1e4)

dimnames(co)=list( levels(all_mono$CELL_LINE), levels(all_mono$COMPOUND) )

require(gplots)
heatmap.2(t(co),trace="none",symbreaks=F, key.title="", margins=c(5,8))

top_two_pcs=irlba(co,2,2)
plot(top_two_pcs$v[,1], top_two_pcs$v[,2], col="white", xlab="PC 1", ylab="PC 2")
text(top_two_pcs$v[,1], top_two_pcs$v[,2], colnames(co))

write.table(co,"pca_imputed_mono.txt",quote=F)
#write.table(co,"pca_imputed_mono.txt",quote=F)

mono_dist=as.matrix(dist(t(co)))
mono_dist=mono_dist/median(mono_dist[upper.tri(mono_dist)])
