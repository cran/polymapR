% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exported_functions.R
\name{finish_linkage_analysis}
\alias{finish_linkage_analysis}
\title{Linkage analysis between all markertypes within LG.}
\usage{
finish_linkage_analysis(
  input_type = "discrete",
  marker_assignment,
  dosage_matrix,
  probgeno_df,
  chk,
  marker_combinations = NULL,
  target_parent = "P1",
  other_parent = "P2",
  convert_palindrome_markers = TRUE,
  ploidy,
  ploidy2 = NULL,
  pairing = c("random", "preferential"),
  prefPars = c(0, 0),
  LG_number,
  verbose = TRUE,
  log = NULL,
  ...
)
}
\arguments{
\item{input_type}{Can be either one of 'discrete' or 'probabilistic'. For the former (default), \code{dosage_matrix} must be supplied,
while for the latter \code{probgeno_df} and \code{chk} must be supplied.}

\item{marker_assignment}{A marker assignment matrix with markernames as rownames and at least containing the column \code{"Assigned_LG"}.}

\item{dosage_matrix}{An integer matrix with markers in rows and individuals in columns.}

\item{probgeno_df}{A data frame as read from the scores file produced by function
\code{saveMarkerModels} of R package \code{fitPoly}, or alternatively, a data frame containing the following columns:
\itemize{
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
Probabilities of dosage score '1' etc. (up to max offspring dosage, e.g. P4 for tetraploid population)
}
\item{maxP}{
Maximum genotype probability identified for a particular individual and marker combination
}
\item{maxgeno}{
Most probable dosage for a particular individual and marker combination
}
\item{geno}{
Most probable dosage for a particular individual and marker combination, if \code{maxP} exceeds a user-defined threshold (e.g. 0.9), otherwise \code{NA}}
}}

\item{chk}{Output list as returned by function \code{\link{checkF1}}. This argument is only needed if probabilistic genotypes are used.}

\item{marker_combinations}{A matrix with four columns specifying marker combinations to calculate linkage.
If NULL all combinations are used for which there are rf functions.
Dosages of markers should be in the same order as specified in the names of rf functions.
E.g. if using 1.0_2.0 and 1.0_3.0 types use: \code{matrix(c(1,0,2,0,1,0,3,0), byrow = TRUE, ncol = 4)}}

\item{target_parent}{Character string specifying target parent.}

\item{other_parent}{Character string specifying other parent.}

\item{convert_palindrome_markers}{Logical. Should markers that behave the same for both parents be converted to a workable format for that parent? E.g.: should 3.1 markers be converted to 1.3?}

\item{ploidy}{Ploidy level of parent 1. If parent 2 has the same ploidy level, then also the ploidy level of parent 2.}

\item{ploidy2}{Integer, by default \code{NULL}. If parental ploidies differ, use this to specify the ploidy of parent 2. Note that in cross-ploidy situations, ploidy2 must be smaller than ploidy.}

\item{pairing}{Type of pairing at meiosis, with options \code{"random"} or \code{"preferential"}.}

\item{prefPars}{The estimates for preferential pairing parameters for parent 1 and 2, in range 0 <= p < 2/3. By default this is c(0,0) (so, no preferential pairing).
See the function \code{\link{test_prefpairing}} and the vignette for more details.}

\item{LG_number}{Number of linkage groups (chromosomes).}

\item{verbose}{Should messages be sent to stdout or log?}

\item{log}{Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.}

\item{\dots}{(Other) arguments passed to \code{\link{linkage}}}
}
\value{
Returns a matrix with marker assignments. Number of linkages of 1.0 markers are artificial.
}
\description{
\code{finish_linkage_analysis} is a wrapper for \code{\link{linkage}}, or in the case of probabilistic genotypes, \code{\link{linkage.gp}}.
The function performs linkage calculations between all markertypes within a linkage group.
}
\examples{
\dontrun{
data("screened_data3", "marker_assignments_P1")
linkages_list_P1<-finish_linkage_analysis(marker_assignment=marker_assignments_P1,
                                          dosage_matrix=screened_data3,
                                          target_parent="P1",
                                          other_parent="P2",
                                          convert_palindrome_markers=FALSE,
                                          ploidy=4,
                                          pairing="random",
                                          LG_number=5)
                                          }
}
