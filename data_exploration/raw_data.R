basedir="../Drug_synergy_data/Raw_Data_csv/ch1_training_combinations/"
files=list.files(path=basedir, pattern="*.csv")

concs=matrix(NA,0,6)
viab=matrix(NA,0,6)
for (fn in files) {
  f=read.csv(paste0(basedir,fn) ,header=F, stringsAsFactors=F)
  cl=strsplit(fn,"\\.")[[1]][3]
  
  # agent 1 monotherapy
  concs=rbind(concs, as.numeric(f[2:7,1]))
  viab=rbind(viab, as.numeric(f[2:7,2]))
  rownames(concs)[nrow(concs)]=f[9,2]
  rownames(viab)[nrow(viab)]=cl
  
  # agent 2 monotherapy
  concs=rbind(concs, as.numeric(f[1,2:7]))
  viab=rbind(viab, as.numeric(f[2,2:7]))
  rownames(concs)[nrow(concs)]=f[10,2]
  rownames(viab)[nrow(viab)]=cl
}

require(reshape2)
require(ggplot2)
pdf("mono_curves.pdf", width=12,height=10)
for (ag in unique(rownames(concs))) {
  reps=rownames(concs)==ag
  co=as.data.frame(concs[reps,])
  vi=as.data.frame(viab[reps,])
  vi$cl=rownames(vi)
  vi$id=1:nrow(vi)
  #rownames(co)=rownames(vi)
  rownames(co)=NULL
  mvi=melt(vi, id.vars = c("cl","id"))
  co$cl=rownames(vi)
  co$id=1:nrow(co)
  mco=melt(co, id.vars = c("cl","id"))
  colnames(mco)=c("cl","id","var","conc")
  all(mvi$Var2==mco$Var2)
  mco$vi=mvi$value
  print( ggplot(mco, aes(conc,vi,col=as.factor(cl),linetype=as.factor(id))) + geom_point() + geom_line()+ scale_x_log10() + theme_bw(base_size = 14) + scale_linetype_manual(values = c(rep("solid", max(mco$id)) ), guide = FALSE) + xlab("concentration") + ylab("viability") + ggtitle(ag) )
}
dev.off()
