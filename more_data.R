
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

pdf("more_data.pdf", width=4, height = 4)
ggplot( scores, aes(setup, score) ) + geom_boxplot( outlier.shape = "none" ) + geom_point( size=3, alpha=.3 ) + theme_bw( base_size = 14 )
dev.off()

load("cached_results/subA_lb2_tissue1_seed1.RData")

ggplot(data.frame(x=names(dist), y=sqrt(o$par$eta_sq_cl)), aes(x,y)) + geom_bar(stat="identity") + ggtitle(paste0("Likelihood: ",format(o$value, digits = 3))) + theme_bw(base_size=16) + ylab("importance") + xlab("") + ylim(0,2.3)
