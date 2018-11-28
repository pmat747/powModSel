\name{liketrace}
\alias{liketrace}
\title{Log-likelihood trace plot}
\description{
This function produces trace plots for the log-likelihood at each or specific temperatures.
}
\usage{
liketrace(x, temp = NULL)
}
\arguments{
  \item{x}{A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.}
  \item{temp}{Position of the temperature to be plotted.}
}
\value{
  It produces a set of log-likelihood traceplots.
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\examples{
data(ligoVirgoSim)
liketrace(ligoVirgoSim, temp = c(1)) # Plotting log-likelihoods for the first temperature
}
