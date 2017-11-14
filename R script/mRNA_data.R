###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
seedNo=14159
set.seed(seedNo)
###*****************************


###*****************************
# REQUIRED LIBRARIES
require("tidyverse")
require("stringr")
require("cowplot")
###*****************************


###*********
# Load relevant data
nameDictionary_RNA_Protein <-  read.csv(file = "../Files/nameDictionary_RNA_Protein.csv")
mRNA_untransformed <-  read.csv(file = "../umut/resDf_mrna_trT_set00_StcAllEx_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noFilter_p1Sf_noNorm.csv")
MyData <-  read.csv(file = "../Files/MyData.csv")
###*********


###*********
# Initial editing of files
names(nameDictionary_RNA_Protein)[names(nameDictionary_RNA_Protein) == 'b.'] <- 'b'
names(mRNA_untransformed)[names(mRNA_untransformed) == 'X'] <- 'mRNA_ID'
###*********


###*********
#join untransformed data to nameDictionary by mRNA_id
mRNA_untransformed %>% 
  dplyr::left_join(nameDictionary_RNA_Protein, by = 'mRNA_ID') %>% 
  subset(select = -c(Protein_id, b)) %>% 
  select(mRNA_ID, gene_name, everything())-> UTmNames

UTmNames %>% subset(select = (-c(mRNA_ID))) -> UTmNames
###*********


###*********
#sum gene counts for each media, filter only ribosomal genes, sum ribosomal gene counts for each media, calculate ribosomal fraction 
UTmNames %>% 
  gather(media, count, MURI_016:MURI_140) %>% 
  group_by(media) %>% mutate(total= sum(count)) %>% 
  filter(str_detect(gene_name,regex('^rp'))) %>% 
  group_by(media) %>% mutate(mRNA = sum(count)) %>% 
  mutate(fraction = mRNA/total)-> UTtotalmRNA

UTtotalmRNA %>% 
  subset(select = -c(gene_name,count, total, mRNA )) %>% 
  distinct(fraction) -> UTmfraction 
###*********


#join Mydata (includes growth conditions and generation times) and fraction by sample
UTmfraction %>% left_join(MyData, by = 'media') %>% na.omit() %>% select(-fraction,fraction)-> UTmready

#average the ribosomal fractions of samples that have the same doubling time
UTmready %>% group_by(doublingTimeMinutes) %>% mutate(fraction_avg = (sum(fraction))/n()) -> UTmready

#unique fraction averages
UTmready %>% distinct(fraction_avg, .keep_all = TRUE) -> UTmfinal

#plot of fractions
UTmready %>% ggplot(aes(x=generations_per_hour, y=fraction, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5)

#plot of fraction averages
UTmfinal %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('mRNA Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5)

#get rid of NaCl and MgSO4 stress 
UTmfinal_alt <- UTmfinal[!grepl("stress", UTmfinal$experiment),]
#get rid of NaCl and MgSO4 stress (non-unique fraction_avg)
UTmfinal_alt2 <- UTmready[!grepl("stress", UTmready$experiment),]

UTmfinal_alt2 %>% ggplot(aes(x=generations_per_hour, y=fraction, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5) 

UTmfinal_alt %>% ggplot(aes(x=generations_per_hour, y=fraction_avg, color=carbonSource)) + xlab('Generations per Hour') + ylab('Ribosomal Fraction') + geom_point() + geom_text(aes(label=media),hjust='inward', vjust=1 , size= 5) 
