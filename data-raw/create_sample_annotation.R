#library(proBatch)
library(tidyverse)
library(reshape2)

data("example_proteome", package = 'proBatch')
data("example_sample_annotation", package = 'proBatch')
sample_annotation <- read.csv("~/Dropbox/batch_effects/data/sample_annotation.csv")

example_sample_annotation1 = sample_annotation %>%
  select(FullRunName, MS_batch.final) %>% merge(example_sample_annotation)

#add instruments
example_sample_annotation1 = example_sample_annotation1 %>%
  mutate(instrument = ifelse(MS_batch.final %in% c('Batch_1','Batch_2'), 'Instr_A',
                             ifelse(MS_batch.final %in% c('Batch_3', 'Batch_4', 'Batch_5'),
                                    'Instr_B', 'Instr_C')))

set.seed(1)
example_sample_annotation1$instrument =
  sample(example_sample_annotation1$imaginary_instrument)

