library(tidyverse)
library(readr)
#library(proBatch)

peptide_annotation = read_csv('data/peptide_annotation.csv')

#load the full proteome
proteome <- read_csv("data/OpenSWATH_proteome.csv")

#extract the QTL proteins
QTLproteins = c("Cat", "Glo1", "Cox7a2l", "Lrpprc", "Pm20d1" , "Dhtkd1" , "Haao" ,
                "Tpmt", "Hdhd3" , "D2hgdh" )

#extract a hanful of proteins with overfitting problem
overfitting_prots = c('Mtch1','Cyfip1', 'Nnt', 'Gclm')

#unite proteins for subsampling
proteins_to_subsample = union(QTLproteins, overfitting_prots)

#subsample the selected proteins and iRTs + Bovine
#extract the spike-ins: iRTs and Bovine:
peptides_sub = peptide_annotation %>%
  filter((Gene %in% proteins_to_subsample) |
           grepl('irt', ProteinName, ignore.case = T) |
           grepl('bovine', Gene, ignore.case = T)) %>%
  pull(Peptide)

#sample random peptides
random_pep = sample_random_peptides(proteome)

#unite peptides for subsampling
peptides_sub = union(peptides_sub, random_pep$peptides)

#filter for peptides selected for subsampling
example_proteome = proteome %>% filter(peptide_group_label %in% peptides_sub)

devtools::use_data(example_proteome)
