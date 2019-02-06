#' @title Stepping-stone sampling
#' @name ss
#' @description This function produces a stepping-stone sampling estimate of the marginal likelihood given a set of log-likelihood values at different temperatures.  It can also be used to produce a generalised stepping-stone sampling estimate (see details below).
#' @usage ss(x, temp = NULL)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.  The temperatures must be sorted in an increasing order.
#' @param temp It indicates the temperatures to be used in the analysis, for instance, c(1,3,K) considers the temperatures at those positions, where \eqn{K} is the number of temperatures.  Note that samples from the prior and posterior must be included in the process, even though the latter are not used by the algorithm.  "NULL" stands for all the temperatures in \code{x}.
#' @details Power posterior methods, among them stepping-stone sampling, rely on a set of samples from different transitional distributions, connecting the prior and the posterior distributions, which is defined by the power posterior density
#'          \deqn{p(\theta) \propto L(x|\theta)^{\beta} \pi(\theta), }
#'          where \eqn{\theta} is the parameter vector, \eqn{0 \le \beta \le 1} is the inverse temperature, \eqn{x} is the data, \eqn{p(\theta)} is the power posterior density, \eqn{L(x|\theta)} is the likelihood function, and \eqn{\pi(\theta)} is the prior density.  See more details about the stepping-stone sampling algorithm in Xie et al. (2011).
#'
#'          \code{ss} can be also used to produce a generalised stepping-stone sampling estimate.  In this case, the power posterior density defines a path between a reference distribution \eqn{g} and the posterior and is given by
#'          \deqn{p(\theta) \propto (L(x|\theta)\pi(\theta))^{\beta} g(\theta)^{1-\beta},}
#'          where  \eqn{0 \le \beta \le 1} is the inverse temperature.  The reference distribution needs to be a reasonable approximation of the posterior.  See more details about the generalised stepping-stone sampling algorithm in Fan et al. (2011).
#' @return It produces a stepping-stone sampling estimate.
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @references Fan, Y., Wu, R., Chen, M.-H., Kuo, L., Lewis, P. O., 2011. Choosing among
#'             partition models in Bayesian phylogenetics. \emph{Molecular Biology and Evolution} 28(1), 523--532.
#'             Xie, W., Lewis, P. O., Fan, Y., Kuo, L., Chen, M.-H., 2011. Improving marginal
#'             likelihood estimation for Bayesian phylogenetic model selection.  \emph{Systematic Biology} \bold{60}(2), 150--160.
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' ss(ligoVirgoSim, temp = NULL)
#' }
#' @export
ss = function(x, temp = NULL){

  if( !is.null(temp)){ # selecting certain temperatures (temps)
    #count = stats::aggregate(logL~invTemp, data = x, FUN = length);
    index = which(diff(x$invTemp)!=0);
    index = cbind(c(1, index + 1), c(index, dim(x)[1]));
    index = index[temp,];
    newX  = unlist(apply(index, 1, function(x)seq(x[1],x[2], by = 1)));
    newX  = x[as.character(newX),];
    #rownames(newX) = NULL; # reseting rownames
    x = newX;
  }

  rownames(x) = NULL; # reseting rownames

  count = stats::aggregate(logL~invTemp, data = x, FUN = length);
  # 'count' could be sorted here, increasing order, if 'x' is not ordered
  diff  = diff(count$invTemp); #=diff(unique(x$temperature));#difference between temperatures
  count = count$logL[-dim(count)[1]]; # Number of elements per temperature/no posterior
  delta = rep(diff, times = count); # replicating diff in original data
  N     = sum(count); # = length(delta)
  logL  = x[1:N, 'logL']; # extracting loglike
  loglDelta = logL * delta; # loglike x delta
  invTemp   = x[1:N, 'invTemp']; # extracting temperature
  Rss   = cbind(invTemp, loglDelta);
  Rss   = stats::aggregate(loglDelta~invTemp, data = Rss, FUN = logplusvec);
  Rss   = Rss["loglDelta"] - log(count); # Individual rate estimates

  return(sum(Rss));

}
