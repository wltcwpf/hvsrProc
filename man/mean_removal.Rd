% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mean_removal}
\alias{mean_removal}
\title{Mean removal}
\usage{
mean_removal(ts, nDC = -1L)
}
\arguments{
\item{ts}{An array of time series}

\item{nDC}{a flag number, if negative, DC is based on the whole windowed record;
if positive, DC is based on the first nDC points.}
}
\value{
The DC time series.
}
\description{
Remove the mean
}
