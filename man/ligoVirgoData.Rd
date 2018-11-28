\name{ligoVirgoData}
\docType{data}
\alias{ligoVirgoData}
\title{LIGO-VIRGO simulations}
\description{
The original dataset contained 4s of simulated Gaussian noise, generated using the design sensitivity curves of two Advanced LIGO detectors (Hanford, Livingston) and the Advanced Virgo detector, plus the GW signal.  \code{ligoVirgoData} contains the likelihood values from samples from power posterior density at different temperatures.
}
\usage{ligoVirgoData}
\format{ A data frame with two variable: \code{logL} which contains log-likelihood values and \code{invTemp} which contains the inverse temperatures.}
\source{<http://www.astro.gla.ac.uk/~jveitch/mcmc/PTMCMC/nest_and_mcmc/>}
\references{
  Maturana Russel, P., Meyer, R., Veitch, J., and Christensen, N. 2018. The stepping-stone algorithm for calculating the evidence of gravitational wave models. arXiv preprint arXiv:1810.04488
}
\keyword{datasets}
