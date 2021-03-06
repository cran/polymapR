% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exported_functions.R
\name{consensus_LG_names}
\alias{consensus_LG_names}
\title{Find consensus linkage group names}
\usage{
consensus_LG_names(
  modify_LG,
  template_SxS,
  modify_SxS,
  merge_LGs = TRUE,
  log = NULL
)
}
\arguments{
\item{modify_LG}{A \code{data.frame} with markernames, linkage group (\code{"LG"}) and homologue (\code{"homologue"}), in which the linkage group numbers will be modified}

\item{template_SxS}{A file with assigned markers of which (at least) part is present in both parents of the template parent.}

\item{modify_SxS}{A file with assigned markers of which (at least) part is present in both parents of the parent of which linkage group number are modified.}

\item{merge_LGs}{Logical, by default \code{TRUE}. If \code{FALSE}, any discrepency in the number of linkage groups will not be merged, but removed instead.
This can be needed if the number of chromosomes identified is not equal between parents, and the user wishes to proceed with a core set.}

\item{log}{Character string specifying the log filename to which standard output should be written. If NULL log is send to stdout.}
}
\value{
A modified modified_LG according to the template_SxS linkage group numbering
}
\description{
Chromosomes that should have same number, might have gotten different numbers between parents during clustering.
\code{consensus_LG_names} uses markers present in both parents (usually 1.1 markers) to modify the linkage group numbers in one parent with the other as template
}
\examples{
data("LGHomDf_P2_2", "P1_SxS_Assigned", "P2_SxS_Assigned")
consensus_LGHomDf<-consensus_LG_names(LGHomDf_P2_2, P1_SxS_Assigned, P2_SxS_Assigned)
}
