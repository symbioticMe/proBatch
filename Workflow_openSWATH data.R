# Workflow with Allruns_requant from OpenSWATH output data 

# OpenSWATH Data 
Allruns_requant = read.table(file = "S:/html/AllRuns_Evosep_feature_alignment_requant.tsv", sep = '\t', header = TRUE)
#Allruns_requant_test = sample_n(Allruns_requant, 10000)
#Allruns_requant_test = Allruns_requant_test %>% mutate(FullRunName = gsub('wevan_(.+)\\.mzXML\\.gz', '\\1', filename))
  # the FUllRunName still doesn't match with example_sample_annotation file in package (order of features differ). 

################### Reference data frames ##################
load("~/R/proBatch/data/example_peptide_annotation.RData")
load("~/R/proBatch/data/example_sample_annotation1.rda")
load("~/R/proBatch/data/example_proteome.RData")
############################################################

############## start from OpenSWATH data ##################
openSWATH_data = Allruns_requant

# Generate df_long 
SWATH_long = openSWATH_data %>% 
  dplyr::select(c("peptide_group_label", "Intensity", "FullRunName"))


# Generate peptide annotation data frame 
peptide_annotation = create_peptide_annotation(openSWATH_data, 
                                               peptide_col = "peptide_group_label", 
                                               annotation_col = c("transition_group_id", "filename" , "RT" ,  "FullPeptideName"  ,  "Charge" ,             
                                                               "m.z" ,  "Intensity" , "ProteinName" ,"decoy", "assay_rt" ,"delta_rt"     ,                 
                                                                "norm_RT"  , "rt_score"  ,  "sn_ratio"  ,  "xx_lda_prelim_score"  ,
                                                               "xx_swath_prelim_score" ,         "aggr_Peak_Area" ,    "aggr_Peak_Apex"    ,            
                                                              "aggr_Fragment_Annotation",       "peak_group_rank"          ,      "d_score"        ,               
                                                              "transition_group_id_m_score"  ,  "transition_group_id_p_value" ,"ProteinName_m_score"     ,      
                                                              "ProteinName_p_value"  ,          "align_runid"   ,                 "align_origfilename"  ,  "align_clusterid"  )) 


# retrieve sample annotation data frame + add order 
sample_annotation = example_sample_annotation1 # PLEASE ADD ONE WITH MATCHING SAMPLE ANNOTATION 
sample_annotation = date_to_sample_order (sample_annotation,
                                                time_column = c('RunDate','RunTime'),
                                                new_time_column = 'DateTime',
                                                dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                                order_col = 'order',
                                                instrument_col = NULL)


############## data pre-processing scripts ##################
# remove peptides with missing batch 
SWATH_long = remove_peptides_with_missing_batch(SWATH_long, sample_annotation,
                                                            batch_col = 'MS_batch.final',
                                                            feature_id_col = 'peptide_group_label')

# clean_requants - requant not necessary in openSWATH file 
SWATH_long = clean_requants(SWATH_long, sample_annotation, peptide_annotation,
                                          batch_col = 'MS_batch.final',
                                          feature_id_col = 'peptide_group_label',
                                          m_score = "transition_group_id_m_score" ,
                                          missing_frac_batch = .3, missing_frac_total = .3)
  #not sure if transition_group_id_m_score is the right column to target!!

############## normalization ###############################
# Generate data_matrix from SWATH_long
SWATH_matrix = convert_to_matrix(SWATH_long, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')


# log2 transformation of data_matrix 
SWATH_matrix_log2 = log2(SWATH_matrix + 1) 

# quantile normalization of the log2 transformed data 
SWATH_matrix_qnorm = quantile_normalize(data_matrix_log2)

# Normalize with LOESS fitting 
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
  #does sample_annotation and samples in matrix need to be in the same order? 

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



