names(nameDictionary_RNA_Protein)[names(nameDictionary_RNA_Protein) == 'b#'] <- 'b'

names(resDf_mrna_trT_set00_StcAllEx_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst)[names(resDf_mrna_trT_set00_StcAllEx_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst) == 'X1'] <- 'mRNA_ID'

resDf_mrna_trT_set00_StcAllEx_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst %>% left_join(nameDictionary_RNA_Protein, by = 'mRNA_ID') %>% subset(select = -c(gene_id, b)) %>% select(mRNA_ID, gene_name, everything())-> mRNA

mRNA %>% subset(select = (-c(mRNA_ID))) -> mRNA

mRNA %>% gather(media, count, MURI_016:MURI_140) %>% group_by(media) %>% mutate(total= sum(count)) %>% filter(str_detect(gene_name,regex('^rp'))) %>% group_by(media) %>% mutate(mrna = sum(count)) %>% mutate(fraction = mrna/total)-> total_mrna

total_mrna %>% subset(select = -c(gene_name,count, total, mrna )) %>% distinct(fraction) -> mrna_fraction 

metaProtein %>% select(growthPhase,carbonSource, dataSet,growthTime_hr,doublingTimeMinutes,Mg_mM_Levels,Na_mM_Levels) %>% filter(growthPhase == 'exponential') -> df2 

df2 %>% mutate(generations_per_hour = ((growthTime_hr*60)/doublingTimeMinutes)/growthTime_hr) -> df2

write.csv(df2, file = "MyData.csv")

mrna_fraction %>% left_join(MyData, by = 'media') %>% na.omit() %>% select(-fraction,fraction)-> ready

ready %>% group_by(doublingTimeMinutes) %>% mutate(fraction_avg = (sum(fraction))/n()) -> ready

ready %>% distinct(fraction_avg, .keep_all = TRUE) -> mrna_final

mrna_final %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource, shape=Na_mM_Levels)) + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5)


