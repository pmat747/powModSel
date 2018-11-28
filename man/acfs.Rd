\name{acfs}
\alias{acfs}
\title{Autocorrelation function}
\description{
The function \code{acfs} computes estimates of the autocorrelation function via \code{acf} function.
}
\usage{
acfs(x, temp = NULL, lag.max = NULL)
}
\arguments{
  \item{x}{A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.}
  \item{temp}{Position of the temperature to be plotted.}
  \item{lag.max}{maximum lag at which to calculate the autocorralation function.}
}
\value{
  It produces a set of autocorrelation plots.
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\examples{
data(ligoVirgoSim)
acfs(ligoVirgoSim, temp = 1:3)
}
