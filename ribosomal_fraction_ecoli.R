library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(stringr)

names(nameDictionary_RNA_Protein)[names(nameDictionary_RNA_Protein) == 'b#'] <- 'b'
names(resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst)[names(resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst) == 'X1'] <- 'gene_id'

#join normalized data to nameDictionary by gene_id
resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst %>% left_join(nameDictionary_RNA_Protein, by = 'gene_id') %>% subset(select = -c(mRNA_ID, b)) %>% select(gene_id, gene_name, everything())-> df

df %>% subset(select = (-c(gene_id))) -> df

#sum gene counts for each media, filter only ribosomal genes, sum ribosomal gene counts for each media, calculate ribosomal fraction 
df %>% gather(media, count, MURI_016:MURI_140) %>% group_by(media) %>% mutate(total= sum(count)) %>% filter(str_detect(gene_name,regex('^rp'))) %>% group_by(media) %>% mutate(ribo = sum(count)) %>% mutate(fraction = ribo/total)-> total_ribo

total_ribo %>% subset(select = -c(gene_name,count, total, ribo )) %>% distinct(fraction) -> fraction 

#select growth condidtions and doubling times, and filter only exponential phase samples
metaProtein %>% select(growthPhase,carbonSource, dataSet,growthTime_hr,doublingTimeMinutes,Mg_mM_Levels,Na_mM_Levels) %>% filter(growthPhase == 'exponential') -> df2 

#calculate generations per hour
df2 %>% mutate(generations_per_hour = ((growthTime_hr*60)/doublingTimeMinutes)/growthTime_hr) -> df2

write.csv(df2, file = "MyData.csv")

#join Mydata (includes growth conditions and generation times) and fraction by sample
fraction %>% left_join(MyData, by = 'media') %>% na.omit() %>% select(-fraction,fraction)-> graphready

#average the ribosomal fractions of samples that have the same doubling time
graphready %>% group_by(doublingTimeMinutes) %>% mutate(fraction_avg = (sum(fraction))/n()) -> graphready

graphready %>% distinct(fraction_avg, .keep_all = TRUE) -> final

final %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5)
 
final_alt <- final[!grepl("stress", final$experiment),]

final_alt %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5) 




