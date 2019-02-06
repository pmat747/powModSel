#' @title Autocorrelation function
#' @description This function produces the autocorrelation plot for the log-likelihood at each temperature
#' @name acfs
#' @usage acfs(x, temp = NULL, lag.max = NULL)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.
#' @param temp Position of the temperature to be plotted.
#' @param lag.max maximum lag at which to calculate the autocorralation function.
#' @return It produces a set of autocorrelation plots.
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' acfs(ligoVirgoSim, temp = 1:3)
#' }
#' @export
acfs = function(x, temp = NULL, lag.max = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));
  j     = dim(index)[1];

  if(is.null(temp)){

    j    = dim(index)[1];
    temp = 1:j;
  }

  tempvalue = round(unique(x$invTemp), 4);

  graphics::par(mar = rep(4,4));

  for(i in temp){
    stats::acf(x$logL[index[i,1]:index[i,2]],
               main = paste( "Inverse Temperature ", tempvalue[i], " - ", i),
               lag.max = lag.max);
  }

  cat("Number of autocorrelation plots", length(temp));
}
