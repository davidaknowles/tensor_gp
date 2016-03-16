
mean(train$SYNERGY_SCORE[train$COMPOUND_A %in% c("AKT","AKT_1") & train$COMPOUND_B %in% c("AKT","AKT_1") ])
hist(train$SYNERGY_SCORE[train$COMPOUND_A %in% c("AKT","AKT_1") & train$COMPOUND_B %in% c("AKT","AKT_1") ])

ss=train$SYNERGY_SCORE[train$COMPOUND_A %in% c("AKT","AKT_1") & train$COMPOUND_B %in% c("AKT","AKT_1") ]
names(ss)=train$CELL_LINE[train$COMPOUND_A %in% c("AKT","AKT_1") & train$COMPOUND_B %in% c("AKT","AKT_1")]
qplot(names(ss),ss)+theme_bw(base_size = 18) +theme(axis.text.x = element_text(angle = 90, hjust = 1))

ALK_IGFR