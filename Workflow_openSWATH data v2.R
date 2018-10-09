# Workflow with Allruns_requant from OpenSWATH output data 

# OpenSWATH Data 
Allruns_requant = read.table(file = "S:/html/AllRuns_Evosep_feature_alignment_requant.tsv", sep = '\t', header = TRUE)
Allruns_requant = Allruns_requant %>% mutate(FullRunName = gsub('wevan_(.+)\\.mzXML\\.gz', '\\1', filename))

peptide_annotation <- read.csv("C:/Users/Chloe/Desktop/Aebersold Lab/2. Batch correction/Data/peptide_annotations_6600evosep.csv")
sample_annotation <- read.csv("C:/Users/Chloe/Desktop/Aebersold Lab/2. Batch correction/Data/sample_annotation_6600evosep.csv")


Fullruns_requant = read.table(file = "S:/html/E1801171630_feature_alignment_requant.tsv", sep = '\t', header = TRUE)
Fullruns_requant = Fullruns_requant %>% mutate(FullRunName = gsub('wevan_(.+)\\.mzXML\\.gz', '\\1', filename))

############## pre-processing ##################
# retrieve sample annotation data frame + add order 
sample_annotation = date_to_sample_order (sample_annotation,
                                                time_column = c('RunDate','RunTime'),
                                                new_time_column = 'DateTime',
                                                dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                                order_col = 'order',
                                                instrument_col = NULL)

# remove peptides with missing batch 
Fullruns_requant = remove_peptides_with_missing_batch(Fullruns_requant, sample_annotation,
                                                            batch_col = 'MS_batch.final',
                                                            feature_id_col = 'peptide_group_label')



############## normalization ###############################
openSWATH_data = Fullruns_requant

# Generate data_matrix from SWATH_long
SWATH_matrix = convert_to_matrix(openSWATH_data, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')


# log2 transformation of data_matrix 
SWATH_matrix_log2 = log2(SWATH_matrix + 1) 

# quantile normalization of the log2 transformed data 
SWATH_matrix_qnorm = quantile_normalize(SWATH_matrix_log2)

# LOESS fitting 
SWATH_list_fit = normalize_custom_fit(SWATH_matrix_qnorm, sample_annotation,
                                      batch_col = 'MS_batch.final',
                                      feature_id_col = 'peptide_group_label',
                                      sample_id_col = 'FullRunName',
                                      measure_col = 'Intensity',
                                      sample_order_col = 'order',
                                      fit_func = fit_nonlinear,
                                      fitFunc = 'loess_regression')
SWATH_matrix_fit = SWATH_matrix_fit$data_matrix

# ComBat noramlization 
SWATH_matrix_ComBat = correct_with_ComBat (SWATH_matrix_fit, sample_annotation,
                                                      batch_col = 'MS_batch.final', par.prior = TRUE)



################### convert to data_long for plotting #########################
# matrix: SWATH_matrix_log2, SWATH_matrix_qnorm, SWATH_matrix_fit, SWATH_matrix_ComBat

SWATH_long_log2 = matrix_to_long(SWATH_matrix_log2, feature_id_col = 'peptide_group_label',
                            measure_col = 'Intensity', sample_id_col = 'FullRunName')

SWATH_long_qnorm = matrix_to_long(SWATH_matrix_qnorm, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity', sample_id_col = 'FullRunName')

SWATH_long_fit = matrix_to_long(SWATH_matrix_fit, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity', sample_id_col = 'FullRunName')

SWATH_long_ComBat = matrix_to_long(SWATH_matrix_ComBat, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity', sample_id_col = 'FullRunName')


################## plotting diagnostics of normalization ###########################



