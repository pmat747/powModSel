\name{ligoVirgoSim}
\docType{data}
\alias{ligoVirgoSim}
\title{LIGO-VIRGO simulations}
\description{
The original dataset contained 4s of simulated Gaussian noise, generated using the design sensitivity curves of two Advanced LIGO detectors (Hanford, Livingston) and the Advanced Virgo detector, plus a gravitational wave signal.  \code{ligoVirgoSim} contains log-likelihood values from samples drawn from a power posterior density at different temperatures (see more details in Maturana R. et al. (2018)).
}
\usage{ligoVirgoSim}
\format{ A data frame with two variable: \code{logL} which contains log-likelihood values and \code{invTemp} which contains the inverse temperatures.}
\source{http://www.astro.gla.ac.uk/~jveitch/mcmc/PTMCMC/nest_and_mcmc/}
\references{
  Maturana Russel, P., Meyer, R., Veitch, J., and Christensen, N. 2018. The stepping-stone algorithm for calculating the evidence of gravitational wave models. arXiv preprint arXiv:1810.04488
}
\keyword{datasets}
