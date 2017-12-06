#required library

library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(readr)
library(stringr)
library(readxl)

Quantitative_Proteomics <- read_excel("~/Desktop/Quantitative Proteomics.xlsx")

Gene_Sectors <- read_excel("~/Desktop/Gene Sectors.xlsx")


#*******************
Gene_Sectors %>% 
  gather(Sector, Gene, `C`:`O`) %>% na.omit() -> Sector

#*******************
#joined Sector and qp df to include only genes listed on supp. material 4
Sector %>% left_join(qp, by= 'Gene') %>% 
  subset(select = -c(Sector.x)) %>% 
  na.omit() -> gs

#*******************
#Made Count into numeric values and created a column with total proteome values
gs$Count <- as.numeric(as.character(gs$Count))

gs %>% group_by(Growth_Rate) %>% 
  mutate(proteome = sum(Count)) -> gs

#*******************
#Filtered only genes that pertain to ribosomal structure and created a ribosomal fraction
gs %>% filter(str_detect(Gene,regex('^rp'))) %>% 
  group_by(Growth_Rate) %>% 
  mutate(ribosomal = sum(Count)) %>% 
  mutate(ribosomal_fraction = ribosomal/proteome)-> gs_fraction

#*******************
#Created a dataframe that includes only the ribosomal fraction and growth rate values
gs_fraction %>% 
  subset(select = -c(Gene,Count, proteome, ribosomal )) %>% 
  distinct(ribosomal_fraction) -> gsfraction 

gsfraction %>% ggplot(aes(Growth_Rate, ribosomal_fraction)) + geom_point()
