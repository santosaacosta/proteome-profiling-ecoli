#required library

library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(readr)
library(stringr)
library(readxl)

Quantitative_Proteomics <- read_excel("~/Desktop/Quantitative Proteomics.xlsx")

#*******************

Quantitative_Proteomics %>% 
  gather(Growth_Rate, Count, `0.45205`:`1.0397`) -> qp

#*******************
#Remove rows without numeric values for genes
qp %>% 
  group_by(Growth_Rate) %>% 
  filter(!(Count == NaN)) -> qp

#*******************
#Made Count into numeric values and created a column with total proteome values
qp$Count <- as.numeric(as.character(qp$Count))

qp %>% group_by(Growth_Rate) %>% 
  mutate(proteome = sum(Count)) -> qp

#*******************
#Filtered only genes that pertain to ribosomal structure and created a ribosomal fraction
qp %>% filter(str_detect(Gene,regex('^rp'))) %>% 
  group_by(Growth_Rate) %>% 
  mutate(ribosomal = sum(Count)) %>% 
  mutate(ribosomal_fraction = ribosomal/proteome)-> qp_fraction

#*******************
#Created a dataframe that includes only the ribosomal fraction and growth rate values
qp_fraction %>% 
  subset(select = -c(Gene,Count, proteome, ribosomal )) %>% 
  distinct(ribosomal_fraction) -> fraction 

fraction %>% ggplot(aes(Growth_Rate, ribosomal_fraction)) + geom_point()
