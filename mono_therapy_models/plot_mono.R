source("mono_therapy_models/mono_model.R")

dc=dcast(mono, CELL_LINE~COMPOUND, value.var = "Einf")
rownames(dc)=dc$CELL_LINE
dc$CELL_LINE=NULL
heatmap(dc)

monoMean=mono %>% group_by(COMPOUND,CELL_LINE) %>% summarize(meanEinf=mean(Einf), replicates=length(Einf))

ggplot(monoMean, aes(COMPOUND,CELL_LINE,fill=meanEinf)) + geom_tile() +  scale_fill_gradient2()  + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(monoMean, aes(COMPOUND,CELL_LINE,fill=replicates)) + geom_tile() +  scale_fill_gradient2()  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
