#### install dependencies ####
#library(ProjectTemplate)
#load.project()

source("http://bioconductor.org/biocLite.R") # from README.md
bioc_deps <- c("GO.db", "preprocessCore", "impute", "sva", "pvca")
cran_deps <- c("tidyverse", "reshape2", "lazyeval", "readr", "WGCNA", "rlang", "corrplot", "ggfortify", "pheatmap")
biocLite(bioc_deps); 
install.packages(cran_deps); 
install.packages("roxygen2"); 
install.packages("lubridate"); 

lapply(bioc_deps, require, character.only = TRUE)
lapply(cran_deps, require, character.only = TRUE)
require(roxygen2)
require(lubridate)


#### Load data 
load("~/R/proBatch/data/example_peptide_annotation.RData")
load("~/R/proBatch/data/example_sample_annotation1.rda")
load("~/R/proBatch/data/example_proteome.RData")

#### Pre-processing of the loaded data ####
# need to pre-process the exmaple_proteome by removing the rows with missing values 
# need to combine the matrix with the peptide annotaiton to be able to align the gene names with the 
  # peptide annotations typically the join_data_matrices function is used to combine the matrices for 
  # comparison of before and after normalization

# generate peptide annotation file from the proteome data 
peptide_annotation = create_peptide_annotation(example_proteome, 
                                               peptide_col = "peptide_group_label", 
                                               protein_col = c("Uniprot_ID", "Gene")) 
  ##comment: May need ot input the expmale protein_col dataset with uniprot_ID and Gene names attached 
  ##so that users can know the format of peptide annotation file to input in order to generate the 
  ##"example_peptide_annotation" dataframe. 

# convert proteome into matrix - works 
data_matrix = convert_to_matrix(example_proteome, feature_id_col = 'peptide_group_label',
                  measure_col = 'Intensity',
                  sample_id_col = 'FullRunName')

# convert the matrix into long data(if necessary) - works 
proteome_long = matrix_to_long(data_matrix, feature_id_col = 'peptide_group_label',
               measure_col = 'Intensity', sample_id_col = 'FullRunName')

# join long matrices (if necessary for comparison before and after normalization)
matrices = join_data_matrices (matrices_list = c(data_matrix, data_matrix_qnorm), step = NULL,
                    sample_annotation = example_sample_annotation1, feature_id_col = 'peptide_group_label',
                    measure_col = 'Intensity',sample_id_col = 'FullRunName')

  ## Comment: for simplicity in the downstream analysis, the example_proteome file (df_long) has been merged with
  ## example_peptide_annotation file that contains more information about the peptide and sample annotation files.  
  ## For some of the functions, separate sample/ peptide annotation files are required while for some, information
  ## are expected to be contained in the df_long (need to double check). 

example_proteome_wPeptide = merge(example_proteome, example_peptide_annotation,
                                  by.x = "peptide_group_label" , by.y = "peptide_group_label")
example_proteome_wSample = merge(example_proteome, example_sample_annotation1,
                                          by.x = "FullRunName" , by.y = "FullRunName")
example_proteome_wPeptide_wSample = merge(example_proteome_wPeptide, example_sample_annotation1,
                                 by.x = "FullRunName" , by.y = "FullRunName")

  ## for the following analysis with pre-processing and analysis, the example_proteome_wPeptide_wSample utilized. 
  ## especially for the functions that require annotations within input dataframe 
  ## (e.g. remove_peptides_with_missing_batch)

# remove peptides with missing batch - works 
example_proteome_wPeptide_wSample = remove_peptides_with_missing_batch(example_proteome_wPeptide_wSample)

# clean requants - works 
example_proteom_requants = clean_requants(example_proteome, example_sample_annotation1, 
                           batch_col = 'MS_batch.final',
                           feature_id_col = 'peptide_group_label',
                           missing_frac_batch = .3, missing_frac_total = .3)

# make order in sample_annotation file - works (fixed)
sample_annotation_order = date_to_sample_order (sample_annotation = example_sample_annotation1,
                                                time_column = c('RunDate','RunTime'),
                                                new_time_column = 'DateTime',
                                                dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                                order_col = 'order',
                                                instrument_col = NULL)
  
  ## Error in mutate_impl(.data, dots) : 
  ## Evaluation error: 'as.POSIXct' is not an exported object from 'namespace:lubridate'. 
  dates_to_posix <- function(sample_annotation,
                             time_column = c('RunDate','RunTime'),
                             new_time_column = NULL,
                             dateTimeFormat = c("%b_%d", "%H:%M:%S")){
    if (length(time_column) == 1){
      if(is.null(new_time_column)) new_time_column = time_column
      time_col = as.character(sample_annotation[[time_column]])
      sample_annotation[[new_time_column]] = as.POSIXct(time_col ,
                                                        format=dateTimeFormat)
    } else {
      sample_annotation = sample_annotation %>%
        mutate(dateTime = paste(!!!syms(time_column), sep=" ")) %>%
        #mutate(dateTime = lubridate::as.POSIXct(dateTime,
        #                                        format = paste(dateTimeFormat, collapse = ' '))) %>%
        mutate(dateTime = as.POSIXct(dateTime, format = paste(dateTimeFormat, collapse = ' '))) %>%        
        rename(!!new_time_column := dateTime)
    }
    return(sample_annotation)
  }
  
  

#### Data transformation, normalization and batch correction ####
# log2 transformation of data_matrix 
data_matrix_log2 = log2(data_matrix + 1) 

# quantile normalization of the log2 transformed data 
data_matrix_log2qnorm = quantile_normalize(data_matrix_log2)

# Normalize with LOESS fitting 
log2qnormfit = normalize_custom_fit(data_matrix_log2qnorm, sample_annotation_order,
                                                batch_col = 'MS_batch.final',
                                                feature_id_col = 'peptide_group_label',
                                                sample_id_col = 'FullRunName',
                                                measure_col = 'Intensity',
                                                sample_order_col = 'order',
                                                fit_func = fit_nonlinear,
                                                fitFunc = 'loess_regression')
data_matrix_log2qnormfit = log2qnormfit$data_matrix
  ## after the LOESs fitting, rownames (peptides) are lost from the matrix. Confirm if anything done wrong 
  ## rownames were added manually from the matrix it was normalized from. 
  rownames(data_matrix_log2qnormfit) = rownames(data_matrix_log2qnorm)


# Batch correction with LOESS+ ComBat 
data_matrix_log2qnormfitCombat = correct_with_ComBat (data_matrix_log2qnormfit, sample_annotation_order,
                                          batch_col = 'MS_batch.final', par.prior = TRUE)

# Batch correction with ComBat ONLY
data_matrix_log2qnormCombat = correct_with_ComBat (data_matrix_log2qnorm, sample_annotation_order,
                                                   batch_col = 'MS_batch.final', par.prior = TRUE)

#### Plot diagnostics for each of the steps ####
  ## for the following analysis, the plots between log2(1), log2+qnorm(2), log2+qnorm+LOESS+Combat(3) and 
  ## Log2+qnorm+ComBat(4) were compared side by side. For simplicity, there were renamed
  data_1 = data_matrix_log2
  data_2 = data_matrix_log2qnorm
  data_3 = data_matrix_log2qnormfitCombat
  data_4 = data_matrix_log2qnormCombat
  
  ## the matrix has been converted back to df_long for some diagnostics that require long data frames 
  df_long_1 = matrix_to_long(data_matrix_log2, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity', sample_id_col = 'FullRunName' )
  df_long_2 = matrix_to_long(data_matrix_log2qnorm, feature_id_col = 'peptide_group_label',
                            measure_col = 'Intensity', sample_id_col = 'FullRunName')
  df_long_3 = matrix_to_long(data_matrix_log2qnormfitCombat, feature_id_col = 'peptide_group_label',
                             measure_col = 'Intensity', sample_id_col = 'FullRunName')
  df_long_4 = matrix_to_long(data_matrix_log2qnormCombat, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')
  
# plot sample mean - works
  ## suggestion for plots: if the users can change the y-axis of the graph, would be great. 
  ## Maybe helpeful for the comparions (look at qnorm data for reference)
plot_sample_mean(data_1, sample_annotation = sample_annotation_order, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(data_2, sample_annotation = sample_annotation_order, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(data_3, sample_annotation = sample_annotation_order, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(data_4, sample_annotation = sample_annotation_order, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)


# plot boxplots - works
plot_boxplot(df_long_1,  sample_annotation = sample_annotation_order,
         sample_id_col = 'FullRunName',
         measure_col = 'Intensity',
         order_col = "order",
         batch_col = 'MS_batch.final',
         facet_col = NULL,
         color_by_batch = T, color_scheme = 'brewer',
         theme = 'classic',
         plot_title = NULL, order_per_facet = F)

plot_boxplot(df_long_2,  sample_annotation = sample_annotation_order,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(df_long_3,  sample_annotation = sample_annotation_order,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(df_long_4,  sample_annotation = sample_annotation_order,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)


# plot clustering - works
colors_list = sample_annotation_to_colors(sample_annotation_order,  columns_for_plotting = NULL,
                                       sample_id_col = 'FullRunName',
                                       factor_columns = c('MS_batch.final','EarTag', "Strain", "Diet", "Sex"),
                                       not_factor_columns = c("Age_Days", "order"),
                                       rare_categories_to_other = T,
                                       numerics_to_log = F,
                                       numeric_palette_type = 'brewer',
                                       granularity = 10)
color_df = colors_list$color_df
plot_sample_clustering(data_1, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)
plot_sample_clustering(data_2, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)
plot_sample_clustering(data_3, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)
plot_sample_clustering(data_4, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)


# plot heatmap - works 
  ## Note: Adding sample annotation file makes an error: 
  ## Error in seq.int(rx[1L], rx[2L], length.out = nb) : 'from' must be a finite number 

plot_heatmap(data_1, sample_annotation = NULL, fill_the_missing = T,
         cluster_rows = T, cluster_cols = F,
         annotation_color_list = colors_list,
         heatmap_color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
         color_for_missing = 'black',
         filename = NA, plot_title = NA)


# plot PCA - works
plot_PCA(data_1, sample_annotation_order, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(data_2, sample_annotation_order, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(data_3, sample_annotation_order, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(data_4, sample_annotation_order, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 


# plot PVCA - works
pvca1 = plot_pvca(data_1, sample_annotation_order,
         sample_id_col = 'FullRunName',
         feature_id_col = 'peptide_group_label',
         technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
         biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
         fill_the_missing = 0,
         threshold_pca = .6, threshold_var = .01,
         colors_for_bars = NULL,
         theme = 'classic', plot_title = NULL)

pvca2 = plot_pvca(data_2, sample_annotation_order,
          sample_id_col = 'FullRunName',
          feature_id_col = 'peptide_group_label',
          technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
          biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
          fill_the_missing = 0,
          threshold_pca = .6, threshold_var = .01,
          colors_for_bars = NULL,
          theme = 'classic', plot_title = NULL)

pvca3 = plot_pvca(data_3, sample_annotation_order,
          sample_id_col = 'FullRunName',
          feature_id_col = 'peptide_group_label',
          technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
          biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
          fill_the_missing = 0,
          threshold_pca = .6, threshold_var = .01,
          colors_for_bars = NULL,
          theme = 'classic', plot_title = NULL)

pvca4 = plot_pvca(data_4, sample_annotation_order,
          sample_id_col = 'FullRunName',
          feature_id_col = 'peptide_group_label',
          technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
          biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
          fill_the_missing = 0,
          threshold_pca = .6, threshold_var = .01,
          colors_for_bars = NULL,
          theme = 'classic', plot_title = NULL)

# plot peptide levels of one protein - works 
plot_peptides_of_one_protein (proteinName = "Haao",  protein_col = "Gene", df_long_1, 
                              sample_annotation_order,example_peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')


# plot peptide trend - Works
  ## Error: Only strings can be converted to symbols (when order_col = NULL)
rows = which(example_peptide_annotation$Gene == "Haao")
peptide_Haao = example_peptide_annotation$peptide_group_label[rows]

plot_peptide_trend (pep_name = peptide_Haao[1], 
                    df_long = df_long_1, 
                    sample_annotation = sample_annotation_order,
                    order_col = "order",
                    sample_id_col = 'FullRunName',
                    batch_col = 'MS_batch.final',
                    measure_col = 'Intensity',
                    feature_id_col = 'peptide_group_label',
                    geom = c('point', 'line'),
                    color_by_batch = F, facet_by_batch = F,
                    requant = NULL,
                    plot_title = NULL,
                    vline_color ='red',
                    theme = 'classic')



# plot spike-ins - Works (fixed)
plot_spike_ins_trend(df_long = df_long_1,
                     sample_annotation = sample_annotation_order,
                     peptide_annotation = example_peptide_annotation,
                     protein_col = 'Gene',
                     order_col = 'order',
                     spike_ins = "BOVINE_A1ag" ,
                     sample_id_col = 'FullRunName',
                     batch_col = 'MS_batch.final',
                     measure_col = 'Intensity',
                     feature_id_col = 'peptide_group_label',
                     requant = NULL,
                     plot_title = 'Spike-in BOVINE protein peptides')

plot_spike_ins_trend(df_long = df_long_1,
                     sample_annotation = sample_annotation_order,
                     peptide_annotation = example_peptide_annotation,
                     protein_col = 'Gene',
                     order_col = 'order',
                     spike_ins = "BOVINE_FetuinB" ,
                     sample_id_col = 'FullRunName',
                     batch_col = 'MS_batch.final',
                     measure_col = 'Intensity',
                     feature_id_col = 'peptide_group_label',
                     requant = NULL,
                     plot_title = 'Spike-in BOVINE protein peptides')

  ## Error in fix.by(by.x, x) : 'by' must specify a uniquely valid column 
  ## it requires the df_long_1 to have "Gene" column already in the data frame to align with the same data frame 
  ## edited the function to remove the error 
  plot_spike_ins_trend <- function(df_long, sample_annotation,
                                   peptide_annotation = NULL,
                                   protein_col = 'ProteinName',
                                   order_col = 'order',
                                   spike_ins = 'BOVIN',
                                   sample_id_col = 'FullRunName',
                                   batch_col = 'MS_batch',
                                   measure_col = 'Intensity',
                                   feature_id_col = 'peptide_group_label',
                                   requant = NULL,
                                   plot_title = 'Spike-in BOVINE protein peptides', ...){
    if (!is.null(peptide_annotation)){
      df_long = df_long %>%
        #merge(peptide_annotation, by = protein_col)
        merge(peptide_annotation)
    }
    spike_in_peptides = df_long %>%
      filter(grepl(spike_ins, !!sym(protein_col))) %>%
      pull(feature_id_col) %>% as.character() %>% unique()
    gg = plot_peptide_trend(spike_in_peptides, df_long = df_long,
                            sample_annotation = sample_annotation,
                            order_col = order_col,
                            sample_id_col = sample_id_col,
                            batch_col = batch_col, measure_col = measure_col,
                            feature_id_col = feature_id_col,
                            plot_title = plot_title, ...)
    return(gg)
  }
  

#### workflow to generate corrletion-based plots ####
#	Plot_protein_corrplot - works but visually difficult to understand 
plot_protein_corrplot(data_matrix = data_1,
                      protein_name = "Haao",
                      peptide_annotation = example_peptide_annotation,
                      protein_col = 'Gene',
                      peptide_col_name = 'peptide_group_label',
                      flavor = 'corrplot',
                      filename = NULL,
                      width = NA, height = NA, unit = c('cm','in','mm'),
                      plot_title = 'peptide correlation matrix')

plot_protein_corrplot(data_matrix = data_3,
                        protein_name = "Haao",
                        peptide_annotation = example_peptide_annotation,
                        protein_col = 'Gene',
                        peptide_col_name = 'peptide_group_label',
                        flavor = 'corrplot',
                        filename = NULL,
                        width = NA, height = NA, unit = c('cm','in','mm'),
                        plot_title = 'peptide correlation matrix')

#	Plot_samples_corrplot  
plot_samples_corrplot(data_matrix = data_1, 
                      samples_to_plot = colnames(data_1)[1:50] ,
                      flavor = 'corrplot', filename = NULL,
                      width = NA, height = NA, unit = c('cm','in','mm'),
                      plot_title = 'Correlation matrix of samples')
  
  
#	Plot_corr_matrix - where to generate corr_matrix? 
#	Get_sample_corr_distrib 
#	Plot_sample_corr_distribution (data_matrix)
#	Get_peptide_corr_df
#	Plot_prot_corr_distribution
#	Plot_within_prot_corr_distribution
#	Plot_peptide_correlation_distr_one_protein (data_matrix_list)







