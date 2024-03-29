% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exported_functions.R
\name{plot_map}
\alias{plot_map}
\title{Plot linkage maps}
\usage{
plot_map(
  maplist,
  highlight = NULL,
  bg_col = "grey",
  highlight_col = "yellow",
  colname_in_mark = NULL,
  colname_beside_mark = NULL,
  palette_in_mark = colorRampPalette(c("white", "purple")),
  palette_beside_mark = colorRampPalette(c("white", "green")),
  color_by_type = FALSE,
  dosage_matrix = NULL,
  parent1 = "P1",
  parent2 = "P2",
  legend = FALSE,
  ...,
  legend.args = list(x = 1, y = 120)
)
}
\arguments{
\item{maplist}{A list of maps. In the first column marker names and in the second their position.}

\item{highlight}{A list of the same length of maplist with vectors of length 2 that specifies the
limits in cM from and to which the plotted chromosomes should be highlighted.}

\item{bg_col}{The background colour of the map.}

\item{highlight_col}{The color of the highlight. Only used if \code{highlight} is specified.}

\item{colname_in_mark}{Optional. The column name of the value to be plotted as marker color.}

\item{colname_beside_mark}{Optional. The column name of the value to be plotted beside the markers.}

\item{palette_in_mark, palette_beside_mark}{Color palette used to plot values. Only used if colnames of the values are specified.}

\item{color_by_type}{Logical. Should the markers be coloured by type? If TRUE, dosage_matrix should be specified.}

\item{dosage_matrix}{Optional (by default \code{NULL}). Dosage matrix of marker genotypes, input of \code{\link{linkage}}}

\item{parent1}{Character string specifying the first (usually maternal) parentname.}

\item{parent2}{Character string specifying the second (usually paternal) parentname.}

\item{legend}{Logical. Should a legend be drawn?}

\item{\dots}{Arguments passed to \code{\link{plot}}}

\item{legend.args}{Optional extra arguments to pass to \code{\link{legend}}, by default a list with x = 1 and y = 120 (position of the legend). 
Additional arguments should be passed using name = value, i.e. as a named list. Note that arguments \code{lty} (= 1) and \code{lwd} (= 2) have
already been used internally (as well as \code{legend} and \code{col}), so cannot be re-specified without causing an error.}
}
\description{
Makes a simple plot of a list of generated linkage maps
}
\examples{
data("maplist_P1")
plot_map(maplist = maplist_P1, colname_in_mark = "nnfit", bg_col = "white",
         palette_in_mark = colorRampPalette(c("blue", "purple", "red")),
         highlight = list(c(20, 60),
         c(60,80),
         c(20,30),
         c(40,70),
         c(60,80)))
}
