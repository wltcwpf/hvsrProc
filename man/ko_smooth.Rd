% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ko_smooth}
\alias{ko_smooth}
\title{Konno Ohmachi smoothing}
\usage{
ko_smooth(freq, amp, b = 20, rate = 2.5)
}
\arguments{
\item{freq}{An array of the frequency}

\item{amp}{An array of the corresponding FAS}

\item{b}{The coefficient of bandwidth. A smaller value will lead to more smoothing}

\item{rate}{The truncation rate of KO smoothing. A larger value will be more accurate but more time consuming}
}
\value{
The smoothed FAS
}
\description{
Apply Konno Ohmachi smoothing to the Fourier Amplitude Spectra (FAS) (Konno, K., Ohmachi, T. Ground-motion characteristics estimated from spectral ratio between horizontal and vertical components of microtremor. Bulletin of the Seismological Society of America. 1998.88 (1): 228–241.)
}
