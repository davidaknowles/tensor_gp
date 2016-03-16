pathways=read.csv("DREAM CHALLENGE TABLE DRUGS-TARGETS LOUKIA.csv",row.names = 1,check.names = F)

table(colSums(pathways))

pathways=pathways[rowSums(pathways)>0,colSums(pathways)>1]

pathways$compound=do.call(rbind,strsplit(rownames(pathways),"-"))[,2]

overl=pathways$compound %in% colnames(pathways)
pathways$compound[ overl ]=paste0(pathways$compound[ overl ]," ")

m=melt(pathways, id.vars = "compound")
m=m[m$value>0,]

allnodes=union(m$compound,m$variable)
n=length(allnodes)

m$compound=factor(m$compound,allnodes)
m$variable=factor(m$variable,allnodes)

#require(Matrix)
#sm=matrix(sparseMatrix(as.numeric(m$compound),as.numeric(m$variable),x=1,dims=c(n,n)))

el=as.matrix(data.frame(from=as.numeric(m$compound), to=as.numeric(m$variable), value=1))
attr(el,"n")=length(allnodes)
vc=(allnodes %in% unique(m$compound))*2 + 2
gplot(el, label=allnodes, gmode="graph", vertex.col="white", label.col =vc, label.pos=5, label.cex=.5,  edge.col = "gray", interactive=T, boxed.labels=F, vertex.cex=0)
