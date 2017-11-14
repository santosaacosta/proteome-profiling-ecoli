#join normalized data to nameDictionary by gene_id
resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst %>% left_join(nameDictionary_RNA_Protein, by = 'gene_id') %>% subset(select = -c(mRNA_ID, b)) %>% select(gene_id, gene_name, everything())-> df

df %>% subset(select = (-c(gene_id))) -> df

#sum gene counts for each media, filter only ribosomal genes, sum ribosomal gene counts for each media, calculate ribosomal fraction 
df %>% gather(media, count, MURI_016:MURI_140) %>% group_by(media) %>% mutate(total= sum(count)) %>% filter(str_detect(gene_name,regex('^rp'))) %>% group_by(media) %>% mutate(ribo = sum(count)) %>% mutate(fraction = ribo/total)-> total_ribo

total_ribo %>% subset(select = -c(gene_name,count, total, ribo )) %>% distinct(fraction) -> fraction 

#select growth condidtions and doubling times, and filter only exponential phase samples
metaProtein %>% select(growthPhase,carbonSource, dataSet,experiment,growthTime_hr,doublingTimeMinutes,Mg_mM,Na_mM,Mg_mM_Levels,Na_mM_Levels) %>% filter(growthPhase == 'exponential') -> df2

#calculate generations per hour
df2 %>% mutate(generations_per_hour = ((growthTime_hr*60)/doublingTimeMinutes)/growthTime_hr) -> df2

write.csv(df2, file = 'df2.csv')
write.csv(MyData, file = "MyData.csv")

#join Mydata (includes growth conditions and generation times) and fraction by sample
fraction %>% left_join(MyData, by = 'media') %>% na.omit() %>% select(-fraction,fraction)-> Na_ready

#average the ribosomal fractions of samples that have the same doubling time
Na_ready %>% group_by(doublingTimeMinutes) %>% mutate(fraction_avg = (sum(fraction))/n()) -> Na_ready

?group_by
?ungroup

Na_ready %>% filter(experiment == 'NaCl_stress') %>% ungroup() %>% distinct(Na_mM, carbonSource,.keep_all = TRUE) -> Na_final

Na_ready %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=Na_mM),hjust='inward', vjust=1 , size= 5) 
Na_final %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=media)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=Na_mM),hjust='inward', vjust=1 , size= 5) 

coef(lm(fraction_avg ~ generations_per_hour, data = Na_final))

#non-distinct ribosomal fraction NaCl stress
Na_ready %>% filter(experiment == 'NaCl_stress') %>% ungroup() -> a
