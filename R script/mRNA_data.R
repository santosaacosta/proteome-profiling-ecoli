####*********
#changing column names 
names(nameDictionary_RNA_Protein)[names(nameDictionary_RNA_Protein) == 'b#'] <- 'b'
names(mRNA_untransformed)[names(mRNA_untransformed) == 'mRNA_id'] <- 'mRNA_ID'
####*********


####*********
#join untransformed data to nameDictionary by mRNA_id
mRNA_untransformed %>% left_join(nameDictionary_RNA_Protein, by = 'mRNA_ID') %>% subset(select = -c(Protein_id, b)) %>% select(mRNA_ID, gene_name, everything())-> UTmNames

UTmNames %>% subset(select = (-c(mRNA_ID))) -> UTmNames
####*********


####*********
#sum gene counts for each media, filter only ribosomal genes, sum ribosomal gene counts for each media, calculate ribosomal fraction 
UTmNames %>% gather(media, count, MURI_016:MURI_171) %>% group_by(media) %>% mutate(total= sum(count)) %>% filter(str_detect(gene_name,regex('^rp'))) %>% group_by(media) %>% mutate(mRNA = sum(count)) %>% mutate(fraction = mRNA/total)-> UTtotalmRNA

UTtotalmRNA %>% subset(select = -c(gene_name,count, total, mRNA )) %>% distinct(fraction) -> UTmfraction 
####*********


####*********
names(metaData_mrna)[names(metaData_mrna) == 'dataSet'] <- 'media'

#join Mydata (includes growth conditions and generation times) and fraction by sample
UTmfraction %>% left_join(metaData_mrna, by = 'media') %>% 
  select(-fraction,fraction)-> UTmready
####*********


####********
#calculate generations per hour
UTmready %>% mutate(generations_per_hour = ((growthTime_hr*60)/doublingTimeMinutes)/growthTime_hr) -> UTmready
####********


####*********
#average the ribosomal fractions of samples that have the same doubling time
UTmready %>% group_by(doublingTimeMinutes) %>% 
  mutate(fraction_avg = (sum(fraction))/n()) -> UTmready
####*********


####*********
#unique fraction averages
UTmready %>% distinct(fraction_avg, .keep_all = TRUE) -> UTmfinal
###**********


####*********
#plot of fractions
UTmready %>% ggplot(aes(x=generations_per_hour, y=fraction, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5)
####*********


####*********
#plot of fraction averages
UTmfinal %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('mRNA Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5)
####*********


####*********
#get rid of NaCl and MgSO4 stress 
UTmfinal_alt <- UTmfinal[!grepl("stress", UTmfinal$experiment),]
#get rid of NaCl and MgSO4 stress (non-unique fraction_avg)
UTmfinal_alt2 <- UTmready[!grepl("stress", UTmready$experiment),]
####*********

UTmfinal_alt2 %>% ggplot(aes(x=generations_per_hour, y=fraction, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5) 

UTmfinal_alt %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5) 
