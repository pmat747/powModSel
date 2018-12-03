#' @title Thermodynamic integration
#' @description It calculates the thermodynamic integration estimate given the likelihood and the temperature values
#' @param x A dataframe
#' @param actPlot Logic value
#' @param temp Integer vector subset of 1:K, where K:number of temperatures
#' @return A numeric value, and an optional plot
#' @return Thermodynamic integration estimate (numeric value)
#' @export ti
#' @importFrom stats

ti = function(x, actPlot = FALSE, temp = NULL){

  if( !is.null(temp)){ # selecting certain temperatures (temp)
    #count = aggregate(logL~invTemp, data = x, FUN = length);
    index = which(diff(x$invTemp)!=0);
    index = cbind(c(1, index + 1), c(index, dim(x)[1]));
    index = index[temp,];
    newX  = unlist(apply(index, 1, function(x)seq(x[1],x[2], by = 1)));
    newX  = x[as.character(newX),];
    #rownames(newX) = NULL; # reseting rownames
    x = newX;
  }

  rownames(x) = NULL; # reseting rownames

  Rti = stats::aggregate(logL~invTemp, FUN = mean, data = x); # Mean per temperature

  if( actPlot == TRUE){
    plot(Rti$invTemp, Rti$logL, xlab = "Inverse temperature",
         ylab = "Log-likelihood mean")
  }

  return( sum( diff(Rti$invTemp) *
          (Rti$logL[-length(Rti$logL)] + Rti$logL[-1]) )/2 );

}
