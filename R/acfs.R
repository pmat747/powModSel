#' @title Autocorrelation plot
#' @description This function produces the autocorrelation plot for the log-likelihood at each temperature
#' @param x dataframe
#' @param temp vector indicating the temperatures to be plotted
#' @param lag.max maximum lag at which to calculate the autocorralation function.
#' @return A set of autocorrelation plot
#' @export acfs

acfs = function(x, temp = NULL, lag.max = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));
  j     = dim(index)[1];

  if(is.null(temp)){

    j    = dim(index)[1];
    temp = 1:j;
  }

  tempvalue = round(unique(x$invTemp), 4);

  par(mar = rep(4,4));

  for(i in temp){
      acf(x$logL[index[i,1]:index[i,2]],
          main = paste( "Temperature ", tempvalue[i], " - ", i),
          lag.max = lag.max);
  }

  cat("Number of autocorrelation plots", length(temp));
}
