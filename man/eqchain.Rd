\name{eqchain}
\alias{eqchain}
\title{Equal chain lengths}
\description{
This function equates the Markov chain lengths at all the temperatures according to the shortest one.
}
\usage{
eqchain(x)
}
\arguments{
\item{x}{A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.}
}
\value{
It produces a new \code{x} with a equal number of log-likelihood values at each temperature.
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\examples{
data(ligoVirgoSim)
eqchain(ligoVirgoSim)
}
