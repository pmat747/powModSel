\name{ss}
\alias{ss}
\title{Stepping-stone sampling}
\description{
This function produces a stepping-stone sampling estimate of the marginal likelihood given a set of likelihood values at different temperatures.  It can also be used to produce a generalised stepping-stone sampling estimate (see details below).
}
\usage{
ss(x, temp = NULL)
}
\arguments{
  \item{x}{A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the likelihood and inverse temperature values, respectively.  The temperatures must be sorted in an increasing order.}
  \item{temp}{It indicates the temperatures to be used in the analysis, for instance, c(1,3,K-1) considers the temperatures at those positions.  Note that samples from the prior must be included in the process.  "NULL" stands for all the temperatures in \code{x}.}
}
\details{
Power posterior methods, among them stepping-stone sampling, rely on a set of samples from different transitional distributions, connecting the prior and the posterior distributions, which is defined by the power posterior density
\deqn{p(\theta) \propto L(x|\theta)^{\beta} \pi(\theta), }
where \eqn{\theta} is the parameter vector, \eqn{0 \le \beta \le 1} is the inverse temperature, \eqn{x} is the data, \eqn{p(\theta)} is the power posterior density, \eqn{L(x|\theta)} is the likelihood function, and \eqn{\pi(\theta)} is the prior density.  See more details about the stepping-stone sampling algorithm in Xie et al. (2011).

\code{ss} can be also used to produce a generalised stepping-stone sampling estimate.  In this case, the power posterior density defines a path between a reference distribution \eqn{g} and the posterior and is given by
\deqn{p(\theta) \propto (L(x|\theta)\pi(\theta))^{\beta} g(\theta)^{1-\beta},}
where  \eqn{0 \le \beta \le 1} is the inverse temperature.  See more details about the generalised stepping-stone sampling algorithm in Fan et al. (2011).
}
\value{
  It produces a stepping-stone sampling estimate.
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\references{

Fan, Y., Wu, R., Chen, M.-H., Kuo, L., Lewis, P. O., 2011. Choosing among
partition models in Bayesian phylogenetics. \emph{Molecular Biology and Evolution} 28(1), 523--532.

Xie, W., Lewis, P. O., Fan, Y., Kuo, L., Chen, M.-H., 2011. Improving marginal
likelihood estimation for Bayesian phylogenetic model selection.  \emph{Systematic Biology} \bold{60}(2), 150--160.

}
\examples{
data(ligovirgodata)
ss(x, temp = NULL)
}
