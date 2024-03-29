% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exported_functions.R
\name{check_maxP}
\alias{check_maxP}
\title{check your dataset's maxP distribution}
\usage{
check_maxP(probgeno_df)
}
\arguments{
\item{probgeno_df}{A data frame as read from the scores file produced by function
\code{saveMarkerModels} of R package \code{fitPoly}, or alternatively, a data frame containing the following columns:
\describe{
\item{SampleName}{
Name of the sample (individual)
}
\item{MarkerName}{
Name of the marker
}
\item{P0}{
Probabilities of dosage score '0'
}
\item{P1...}{
Probabilities of dosage score '1' etc. (up to max dosage, e.g. P4 for tetraploid population)
}
\item{maxP}{
Maximum genotype probability identified for a particular individual and marker combination
}
\item{maxgeno}{
Most probable dosage for a particular individual and marker combination
}
\item{geno}{
Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}
}
}}
}
\value{
This function does not return any value, is simply a visualisation tool to help assess data quality.
}
\description{
Function to assess the distribution of maximum genotype probabilities (\code{maxP}), if these are available. The function
plots a violin graph showing the distribution of the samples' \code{maxP}.
}
\examples{
data("gp_df")
check_maxP(gp_df)
}
