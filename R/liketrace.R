#' @title Log-likelihood trace plot
#' @name liketrace
#' @description This function produces trace plots for the log-likelihood at each or specific temperatures.
#' @usage liketrace(x, temp = NULL)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.
#' @param temp Position of the temperature to be plotted.
#' @return  It produces a set of log-likelihood traceplots.
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' liketrace(ligoVirgoSim, temp = c(1)) # Plotting log-likelihoods for the first temperature
#' }
#' @export
liketrace = function(x, temp = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));

  if(is.null(temp)){

    j    = dim(index)[1];
    temp = 1:j;
  }

  tempvalue  = round(unique(x$invTemp), 4);

  graphics::par(mar = rep(4,4));

  for(i in temp){

    aux = index[i,1]:index[i,2];
    stats::ts.plot(x$logL[aux],
                   main = paste( "Inverse Temperature ", tempvalue[i], " - ", i),
                   ylab = "log-likelihood", xlab = "Iteration");
  }

  cat("Number of traceplots", length(temp));

}
