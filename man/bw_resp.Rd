% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bw_resp.R
\name{bw_resp}
\alias{bw_resp}
\title{Butterworth filter frequency response}
\usage{
bw_resp(freq, fc, nPole, is_causal)
}
\arguments{
\item{freq}{An array of frequency}

\item{fc}{The corner frequency for low-pass/high-pass}

\item{nPole}{The pole parameter for the order, where the order is 2 * nPole. Positive for low-pass; Negative for high-pass}

\item{is_causal}{A binary to indicator if causal is applied. Apply causal if is_causal = TRUE and apply acausal if is_causal = FALSE}
}
\value{
The frequency response
}
\description{
This function calculates Butterworth filter frequency response (for both causal and acausal)
}
