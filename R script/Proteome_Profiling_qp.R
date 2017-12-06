#change the name of columns on data frames to enable joining of tables

names(nameDictionary_RNA_Protein)[names(nameDictionary_RNA_Protein) == 'b#'] <- 'b'
names(nameDictionary_RNA_Protein)[names(nameDictionary_RNA_Protein) == 'Protein_id'] <- 'Gene'
names(resDf_protein)[names(resDf_protein) == 'X1'] <- 'Gene'


#join normalized data to nameDictionary by gene_id

resDf_protein %>% 
  left_join(nameDictionary_RNA_Protein, by = 'Gene') %>% 
  subset(select = -c(mRNA_ID, b)) %>% 
  select(Gene, gene_name, everything())-> df

#remove Gene column, change name of gene_name to Gene, and join df to QP

df %>% subset(select = -c(Gene)) -> df
names(df)[names(df) == 'gene_name'] <- 'Gene'
df %>%
  right_join(Quantitative_Proteomics, by = 'Gene') -> a
 a %>% 
  subset(select = -c(Sector:`1.0397`)) -> a

#remove genes for which there is no data available
a %>% na.omit() -> b

#sum gene counts for each sample, filter only ribosomal genes, sum ribosomal gene counts for each media, calculate ribosomal fraction 
b %>% 
  gather(sample, count, MURI_016:MURI_140) %>% 
  group_by(sample) %>% mutate(total= sum(count)) %>% 
  filter(str_detect(Gene,regex('^rp'))) %>% group_by(sample) %>% 
  mutate(ribo = sum(count)) %>% mutate(fraction = ribo/total)-> total_ribo

total_ribo %>% 
  subset(select = -c(Gene,count, total, ribo )) %>% 
  distinct(fraction) -> fraction 

#change the name of columns on data frames to enable joining of tables
names(MyData)[names(MyData) == 'exp'] <- 'sample'

#join Mydata (includes growth conditions and generation times) & fraction by sample
fraction %>% left_join(MyData, by = 'sample') %>% 
  na.omit() %>% 
  select(-fraction,fraction)-> graphready

#average the ribosomal fractions of samples that have the same doubling time
graphready %>% 
  group_by(doublingTimeMinutes) %>%
  mutate(fraction_avg = (sum(fraction))/n()) -> graphready

#unique fraction averages
graphready %>% distinct(fraction_avg, .keep_all = TRUE) -> final

#plot of fractions (non-averaged)
graphready %>% 
  ggplot(aes(x=generations_per_hour, y=fraction, color=carbonSource)) + 
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=sample),hjust='inward', vjust=1 , size= 5)

#plot of unique fraction averages 
final %>% 
  ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + 
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=sample),hjust='inward', vjust=1 , size= 5)

#get rid of NaCl and MgSO4 stress (non-unique fraction_avg)
final_alt <- graphready[!grepl("stress", graphready$experiment),]
#get rid of NaCl and MgSO4 stress 
final_alt2 <- final[!grepl("stress", final$experiment),]

#plot of ribosomal fractions without salt stress (non-unique fraction_avg)
final_alt %>% ggplot(aes(x=generations_per_hour, y=fraction, color=carbonSource)) + 
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=sample),hjust='inward', vjust=1 , size= 5) 

#plot of ribosomal fractions without salt stress
final_alt2 %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + 
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=sample),hjust='inward', vjust=1 , size= 5) 

#non-distinct ribosomal fraction NaCl stress
graphready %>% filter(experiment == 'NaCl_stress') %>% ungroup() -> Na_non

Na_non %>% ggplot(aes(x=generations_per_hour, y=fraction, color=sample)) +
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=Na_mM),hjust='inward', vjust=1 , size= 5) 

#keep only samples on NaCl stress with unique Na concentrations (5,100,200,300)
graphready %>% 
  filter(experiment == 'NaCl_stress') %>% 
  ungroup() %>% distinct(Na_mM, carbonSource,.keep_all = TRUE) -> Na_final

Na_final %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=sample)) +
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=Na_mM),hjust='inward', vjust=1 , size= 5) 

#non-distinct ribosomal fraction MgSO4 stress
graphready %>% filter(experiment == 'MgSO4_stress_high') %>% ungroup() -> Mg_non

Mg_non %>% ggplot(aes(x=generations_per_hour, y=fraction, color=sample)) + 
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=Mg_mM),hjust='inward', vjust=1 , size= 5) 

#keep only ond samples on MgSO4 stress with unique Mg concentrations (0.08,0.8,8,200)
graphready %>% 
  filter(experiment == 'MgSO4_stress_high') %>% 
  ungroup() %>% distinct(Mg_mM, carbonSource,.keep_all = TRUE) -> Mg_final

Mg_final %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=sample)) + 
  xlab('Generations per Hour') + ylab('Ribosomal Fraction') + 
  geom_point() + geom_text(aes(label=Mg_mM),hjust='inward', vjust=1 , size= 5) 
