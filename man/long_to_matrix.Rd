% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{long_to_matrix}
\alias{long_to_matrix}
\title{Long to wide data format conversion}
\usage{
long_to_matrix(
  df_long,
  feature_id_col = "peptide_group_label",
  measure_col = "Intensity",
  sample_id_col = "FullRunName",
  qual_col = NULL,
  qual_value = 2
)
}
\arguments{
\item{df_long}{data frame where each row is a single feature in a single
sample. It minimally has a \code{sample_id_col}, a \code{feature_id_col}
and a \code{measure_col}, but usually also an \code{m_score} (in OpenSWATH
output result file). See \code{help("example_proteome")} for more details.}

\item{feature_id_col}{name of the column with feature/gene/peptide/protein
ID used in the long format representation \code{df_long}. In the wide
formatted representation \code{data_matrix} this corresponds to the row
names.}

\item{measure_col}{if \code{df_long} is among the parameters, it is the
column with expression/abundance/intensity; otherwise, it is used
internally for consistency.}

\item{sample_id_col}{name of the column in \code{sample_annotation} table,
where the filenames (colnames of the \code{data_matrix} are found).}

\item{qual_col}{column to color point by certain value denoted
by \code{qual_value}. Design with inferred/requant values in
OpenSWATH output data,
which means argument value has to be set to \code{m_score}.}

\item{qual_value}{value in \code{qual_col} to color. For OpenSWATH data,
this argument value has to be set to \code{2} (this is an \code{m_score}
value for imputed values (requant values).}
}
\value{
\code{data_matrix} (\link{proBatch}) like matrix
(features in rows, samples in columns)
}
\description{
Convert from a long data frame representation to a wide matrix representation
}
\examples{
data("example_proteome", package = "proBatch")
proteome_matrix <- long_to_matrix(example_proteome)

}
\seealso{
Other matrix manipulation functions: 
\code{\link{matrix_to_long}()}
}
\concept{matrix manipulation functions}
