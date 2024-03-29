% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/valueIndex.R
\name{valueIndexXD}
\alias{valueIndexXD}
\title{VALUE atomic and spatial index calculation}
\usage{
valueIndexXD(
  ts = NULL,
  dates = NULL,
  index.code = NULL,
  parallel = FALSE,
  max.ncores = 16,
  ncores = NULL
)
}
\arguments{
\item{ts}{list containing time series for each gridpoint/station}

\item{dates}{A character (or \code{POSIXct}) vector following the format \dQuote{YYYY-MM-DD}
(i.e., \code{format = "\%Y-\%m-\%d"} as in \code{\link{strptime}}). Note that the dates are not required 
by all indices, so the default is \code{NULL}. The function will yield an error if dates are required but not supplied.}

\item{index.code}{Characher of the index code to be computed (use VALUE::show.indices).}

\item{parallel}{Use parallel computing?}

\item{max.ncores}{Maximum number of cores to use}

\item{ncores}{Number of cores to use for parallel computing. Default to ncores-1.}
}
\value{
A grid of the index or a list containing the grid of the index and the 
grid of NA percenatage
}
\description{
VALUE index calculation helper function
}
\author{
M. Iturbide, M.N. Legasa
}
