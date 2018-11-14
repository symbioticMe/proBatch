#' Sample annotation data version 1
#'
#' This is data from Evan's aging study with mock instruments to show how
#' instrument-specific functionality works
#'
#' @format A data frame with 375 rows and 18 variables:
#' \describe{
#'   \item{FullRunName}{name of the file, in most functions used for `sample_id_col`}
#'   \item{MS_batch}{mass-spectrometry batch: 7-level factor of manually annotated batches}
#'   \item{EarTag}{mouse ID, i.e. ID of the biological object}
#'   \item{Strain}{mouse strain ID - biological covariate #1}
#'   \item{Diet}{diet - either `HFD` = `High Fat Diet` or `CD` = `Chow Diet`.
#'   `Mix` stands for mixture of several samples}
#'   \item{Sex}{mice sex - 3-level biological covariate. Possible values - ''}
#'   \item{Age_Days}{mice age at sampling - numeric biological covariate}
#'   \item{RunDate}{mass-spectrometry running date. In combination with `RunTime` used for running order determination}
#'   \item{RunTime}{mass-spectrometry running time. In combination with `RunDate` used for running order determination}
#'   \item{SacrificeDate}{date of mouse sacrifice - technical covariate}
#'   \item{ProteinPrepDate}{date of protein preparation - technical covariate}
#'   ...
#' }
"example_sample_annotation"

#' Example protein data
#'
#' This is data from Evan's aging study with all iRT, spike-in peptides,
#' few random peptides and QTL proteins for biological signal improvement demonstration
#'
#' @format A data frame with 200625 rows and 9 variables:
#' \describe{
#'   \item{peptide_group_label}{peptide ID, which is regular feature level. This column is mostly used as `feature_id_col`}
#'   \item{RT}{retention time. Relevant to identify retention time related bias}
#'   \item{Intensity}{peptide group intensity in given sample. Used in function as `measure_col`}
#'   \item{ProteinName}{Protein group ID, specified as N/UniProtID1|UniProtID2|...,
#'   where N is number of protein peptide group maps to. If 1/UniProtID, then this is proteotypic peptide}
#'   \item{Gene}{conventional gene name of corresponding ProteinName}
#'   \item{assay_rt}{retention time as in DIA library}
#'   \item{m_score}{peptide group identification FDR as determined by pyProphet}
#'   \item{FullRunName}{name of the file, in most functions used for `sample_id_col`}
#'   #'   ...
#' }
#' @source PRIDE ID will be added in future
"example_proteome"

#' Peptide annotation data
#'
#' This is data from Evan's aging study annotated with gene names
#'
#' @format A data frame with 200625 rows and 8 variables:
#' \describe{
#'   \item{Peptide}{peptide group label ID, identical to `peptide_group_label` in `example_proteome`}
#'   \item{Gene}{HUGO gene ID}
#'   \item{ProteinName}{protein group name as specified in `example_proteome`}
#'   ... other parameters determined by 'summarize_peptides' function
#'   }
"example_peptide_annotation"
