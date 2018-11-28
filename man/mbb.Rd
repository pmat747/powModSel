\name{mbb}
\alias{mbb}
\title{Moving Block Bootstrap}
\description{
This function produces a set of marginal likelihood estimates for moving block bootstrap observations.  Its main use is for calculating the standard error associated to the thermodynamic integration and stepping-stone sampling estimates.
}
\usage{
mbb(x, l_mbb, n_bootstrap, temp = NULL)
}
\arguments{
  \item{x}{A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.}
  \item{l_mbb}{Block lenghts.}
  \item{n_bootstrap}{Number of bootstrap observations to be analysed.}
  \item{temp}{It indicates the temperatures to be used in the analysis, for instance, c(1,3,K) considers the temperatures at those positions, where \eqn{K} is the number of temperatures.  In this case, the temperatures must be sorted in an increasing order.  Note that samples from the prior and posterior must be included in the process.  \code{NULL} stands for all the temperatures in \code{x}.}
}
\details{
For a block length equal to 1 (\code{l_mbb}=1) the original bootstrap method for i.i.d. data is recovered.  A block lenght greater than one allows to take into account a potential autocorrelation within the Markov chains.  \code{mbb} is being designed to take also into account potential cross-correlation between the Markov chains due to swaps in parallel tempering sampling.  See more details in Maturana R. et al. (2018)
}
\value{
  It produces a list with the following elements:
  \item{Zs}{Marginal likelihood estimates via \code{ti} and \code{ss}.}
  \item{se}{Standard deviation of the marginal likelihood estimates calculated for the bootstrap observations.}
  \item{res}{Marginal likelihood estimate differences between the ones calculated for the bootstrap observations and the original dataset \code{x}.}
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\references{

Kunsch, H. R. 1989. The Jackknife and the Bootstrap for General Stationary Observations. \emph{The Annals of Statistics} \bold{17}(3), 1217--1241.

Maturana Russel, P., Meyer, R., Veitch, J., and Christensen, N. 2018. The stepping-stone algorithm for calculating the evidence of gravitational wave models. arXiv preprint arXiv:1810.04488

}
\examples{
data(ligoVirgoSim)
R = mbb(ligoVirgoSim, l_mbb = 10, n_bootstrap = 20, temp = NULL)
R$se; # standard error of the marginal likelihood estimates
}
