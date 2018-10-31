# FullDataSet Finalized normalization workflow for Evan's analysis #

library(readr) #this is full Evan's dataset
Fullruns_dataset = read_delim("S:/html/E1801171630_feature_alignment_requant.tsv", "\t", escape_double = FALSE)
Fullruns_dataset = Fullruns_dataset %>% mutate(FullRunName = gsub('/scratch/55808263.tmpdir/wevan_(.+)\\.mzXML\\.gz', '\\1', filename))
  # consider using fread if this changes dataset in any way 

# Peptide annotation  
peptide_annotation = read.csv("S:/html/peptide_annotations_6600evosep.csv")

# Sample annotation 
sample_annotation = read.csv("S:/html/sample_annotation_6600evosep.csv")
sample_annotation = date_to_sample_order (sample_annotation,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          order_col = 'order',
                                          instrument_col = NULL)

  # check if any pre-processing step is required for Evan's dataset e.g. requants/ peptide of missing batch etc 

# Generate data_matrix from SWATH_long
SWATH_matrix = convert_to_matrix(Allruns_requant, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')

############################### Normalization & Batch correction ####################################################
# log2 transformation of data_matrix 
SWATH_matrix_log2 = log2(SWATH_matrix + 1) 

# Quantile normalization of the log2 transformed data 
SWATH_matrix_qnorm = quantile_normalize(SWATH_matrix_log2)

# Batch effect correction 
batch_corrected = correct_batch(data_matrix = SWATH_matrix_qnorm, sample_annotation, fitFunc = 'loess_regression', 
                                discreteFunc = 'ComBat', batch_col = 'MS_batch.final',  
                                feature_id_col = 'peptide_group_label', sample_id_col = 'FullRunName',
                                measure_col = 'Intensity',  sample_order_col = 'order', 
                                loess.span = 0.75, abs.threshold = 5, pct.threshold = 0.30)

############################# Preliminary plots to confirm the workflow #############################################

