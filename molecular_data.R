# --------GENE EXPRESSION -----------
require(data.table)
a=fread("Sanger_molecular_data/gex.csv")
setDF(a)
rownames(a)=a$V1
a$V1=NULL
hist(as.numeric(as.matrix(a)))
s=apply(a,1,sd)
hist(s)

b=fread("Sanger_molecular_data/cell_info.csv")
setDF(b)
rownames(b)=b$Sanger.Name

dir.create("processed_data")

tissue_dist=1-outer(b$Disease.Area, b$Disease.Area, FUN = "==")
dimnames(tissue_dist)=list( rownames(b), rownames(b) )
write.csv( tissue_dist, file="processed_data/tissue_dist.csv", quote=F )

setdiff( b$Sanger.Name, colnames(a)  ) # missing cell lines

ccle=fread("CCLE_Expression_Entrez_2012-09-29.gct")
setDF(ccle)
gn=ccle$Description
ccle$Name=NULL
ccle$Description=NULL
ccle=as.matrix(ccle)
rownames(ccle)=gn

#ccle_cell_lines=as.character(read.table("~/Dropbox/ccle_data/cl_names.txt")$V1)
missingCL=b[ setdiff( b$Sanger.Name, colnames(a)  ) , "CCLE.Name" ]
missingCL %in% colnames(ccle)

sanger=as.matrix(a)
colnames(sanger)=b[ colnames(sanger) , "CCLE.Name" ]
sangerMissing=matrix(NA, nrow(sanger), length(missingCL), dimnames = list(rownames(sanger),missingCL))
commonGenes=intersect( rownames(sanger), rownames(ccle) )
commonCL=intersect(colnames(sanger),colnames(ccle))

preds=foreach (g=commonGenes, .combine = rbind) %dopar% {
  geC=ccle[g,commonCL]
  geS=sanger[g,commonCL]
  l=lm(geS ~ geC)
  predict(l, newdata=data.frame( geC=ccle[g,missingCL] ))
}
rownames(preds)=commonGenes

cors=foreach (g=commonGenes, .combine = c) %dopar% {
  geC=ccle[g,commonCL]
  geS=sanger[g,commonCL]
  cor(geS , geC)
}

sangerMissing[commonGenes,]=preds

sangerAll=cbind(sanger,sangerMissing)
sangerAllImp=impute.knn(sangerAll)$data

d=as.matrix(dist(t(sangerAllImp)))
pdf("gex_hclust.pdf",width=15,height=8); plot( hclust( as.dist( d ) ) ); dev.off(); 

d=as.matrix(dist(t(sangerAllImp)))
write.csv( d, file="processed_data/gex_dist.csv", quote=F )

sangerAllImpScaled=t(scale(t(sangerAllImp)))
d=as.matrix(dist(t(sangerAllImpScaled)))
d=d/median(d[upper.tri(d)])
write.csv( d, file="processed_data/gex_scaled_dist.csv", quote=F )

rownames(b)=make.unique(b$CCLE.Name)
colnames(sangerAllImp)=c(colnames(a),b[missingCL,"Sanger.Name"])
write.csv( sangerAllImp, file="processed_data/gex_imputed_with_ccle.csv", quote=F )

#b$Sanger.Name==colnames(d) # not true, GE for two cell lines is missing

tissues=b[colnames(a),"CCLE.Name"]
dimnames(d)=list(tissues,tissues)

# TODO: impute or find missing GE

#------------------ MUTATIONS ---------------------
a=fread("~/Box Sync/astrazeneca_dream/code/Sanger_molecular_data/mutations.csv")
setDF(a)

pdf("mut_counts.pdf",height=12,width=8)
par(mar=c(5,10,5,2)); barplot(table(a$cell_line_name), horiz = T, las=T ); par(mar=c(5,4,4,2))
dev.off()

a$cell_line_name=factor(a$cell_line_name)
a$Gene.name=factor(a$Gene.name)
require(Matrix)
sm=sparseMatrix(i=as.numeric(a$cell_line_name), j=as.numeric(a$Gene.name), x=rep(1,nrow(a)), dimnames = list(levels(a$cell_line_name),levels(a$Gene.name)))
#d=as.matrix(dist(sm))
smbin=t(as.matrix(sm))>0
d=as.matrix(dist(t(smbin)))
require(prabclus)
jd=jaccard(smbin)
dimnames(jd)=list(rownames(sm),rownames(sm))
write.csv( as.matrix(jd),file="mut_dist.csv", quote=F)
dimnames(jd)=list(tissues,tissues)
pdf("mut_hclust_jacc.pdf",width=15,height=8); plot( hclust( as.dist( jd ) ) ); dev.off(); 

#qplot( d[upper.tri(jd)], jd[upper.tri(jd)], col=same_tissue[upper.tri(jd)], xlab="Euclidean", ylab="Jaccard", alpha=.5) + theme_bw()

#qplot( d[upper.tri(jd)], col=same_tissue[upper.tri(jd)], geom="blank" ) + geom_density(position = "identity")
#qplot( jd[upper.tri(jd)], col=same_tissue[upper.tri(jd)], geom="blank" ) + geom_density(position = "identity")

#jd2=d^2 / (d^2 + xy)

#tissues=b[as.character(levels(a$cell_line_name)),"CCLE.Name"]
#dimnames(d)=list(tissues,tissues)
#pdf("mut_hclust.pdf",width=15,height=8); plot( hclust( as.dist( d ) ) ); dev.off(); 

#----------------- CNV ---------------------
a=fread("Sanger_molecular_data/cnv/cnv_gene.csv")
setDF(a)
#qplot(a$min_cn,a$max_cn) + stat_binhex()
require(reshape2)
a$mean_cn=.5*(a$max_cn+a$min_cn)
temp=a[,c("gene","cell_line_name","mean_cn")]
mat=dcast(temp, gene ~ cell_line_name)
rownames(mat)=mat$gene
mat$gene=NULL
mat=as.matrix(mat)
ta=table(mat)
barplot(ta, ylab="count", xlab="copy number")
require(impute)
ik=impute.knn(mat)
matImp=ik$data
cnvdist=as.matrix(dist(t(matImp)))
write.csv( as.matrix(cnvdist),file="cnv_dist.csv", quote=F)

tissues=b[rownames(cnvdist),"CCLE.Name"]
dimnames(cnvdist)=list(tissues,tissues)
pdf("cnv_hclust.pdf",width=15,height=8); plot( hclust( as.dist( cnvdist ) ) ); dev.off(); 

#----------------- MEHTYL ---------------------
a=fread("Sanger_molecular_data/methyl/CpG_probe_level/methyl_probe_beta.csv")
setDF(a)
rownames(a)=a$V1
a$V1=NULL
a=as.matrix(a)

missingCL=setdiff(b$Sanger.Name,colnames(a))

ge=fread("processed_data/gex_imputed_with_ccle.csv")
setDF(ge)
rownames(ge)=ge$V1
ge$V1=NULL
ge=as.matrix(ge)

#temp=rbind( ge[, c(colnames(a),missingCL)], cbind(a, matrix(NA,nrow(a),length(missingCL))))
#temp2=impute.knn(temp, colmax = .99)$data
#methimp=temp2[(nrow(ge)+1):nrow(temp2),]

gedist=as.matrix(dist(t(ge)))
imp=do.call(cbind,foreach (i = missingCL) %do% {
  gedist[i,missingCL]=NA
  o=colnames(gedist)[order(gedist[i,],decreasing = F)[1:10]]
  rowMeans(a[,o])
})
colnames(imp)=missingCL
methimp=cbind(a,imp)

methdist=as.matrix(dist(t(methimp)))
write.csv( as.matrix(methdist),file="processed_data/methyl_dist.csv", quote=F)
#methdist=read.csv("../code/processed_data/methyl_dist.csv",row.names=1,check.names = F)

rownames(b)=b$Sanger.Name
tissues=b[rownames(methdist),"CCLE.Name"]
dimnames(methdist)=list(tissues,tissues)
pdf("meth_hclust.pdf",width=15,height=8); plot( hclust( as.dist( methdist ) ) ); dev.off(); 
