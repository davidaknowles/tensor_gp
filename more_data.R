
sub_challenge="A"
use_tissue=T

test=load_data( "ch1_train_combination_and_monoTherapy.csv","ch1_LB.csv")$test

source("load_cell_line_data.R")
source("load_response_data.R")

cls=levels(train$CELL_LINE)

scores=foreach(setup=c("lb","lb2"), .combine = rbind ) %do% {
    data.frame(setup=setup, score=foreach(i=1:10, .combine=c) %do% {
        resfile=paste0("cached_results/sub",sub_challenge,"_",setup,"_tissue",as.numeric(use_tissue),"_seed",i,".RData")
        if (!file.exists(resfile)) return(NULL)
        load(resfile)
        get_score(o$par$ytest, test)
    } )
}

pdf("more_data.pdf")
ggplot( scores, aes(setup, score) ) + geom_boxplot( outlier.shape = "none" ) + geom_point( size=3, alpha=.3 ) + theme_bw( base_size = 14 )
dev.off()

