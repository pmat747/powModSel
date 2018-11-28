#' @title Log-likelihood trace plot
#' @description This function produces the traceplot for the log-likelihood at each temperature
#' @param x dataframe
#' @param temp vector indicating the temperatures to be plotted
#' @return A set of traceplots
#' @export liketrace

liketrace = function(x, temp = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));

  if(is.null(temp)){

    j    = dim(index)[1];
    temp = 1:j;
  }

  tempvalue  = round(unique(x$invTemp), 4);

  par(mar = rep(4,4));

  for(i in temp){

    aux = index[i,1]:index[i,2];
    ts.plot(x$logL[aux],
            main = paste( "Inverse Temperature ", tempvalue[i], " - ", i),
            ylab = "log-likelihood", xlab = "Iteration");
  }

  cat("Number of traceplots", length(temp));

}
