\name{Burnin}
\alias{burnin}
\title{Burnin period}
\description{
This function discards the first log-likelihood values at each inverse temperature as burnin period.
}
\usage{
burnin(x, percentage = NULL, n_burnin = NULL)
}
\arguments{
  \item{x}{A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.}
  \item{percentage}{Percentage of log-likelihod values to discard at each temperature.}
  \item{n_burnin}{Number of log-likelihod values to discard at each temperature.}
}
\value{
  It produces a new \code{x} with a burnin period.
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\examples{
data(ligoVirgoSim)
burnin(ligoVirgoSim, percentage=10); #discarding 10\% of log-likelihood values at each temperature
}
