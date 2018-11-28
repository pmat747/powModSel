#' @title Log-likelihood trace plot
#' @description This function produces the traceplot for the log-likelihood at each temperature
#' @param x dataframe
#' @param temp vector indicating the temperatures to be plotted
#' @return A set of traceplots
#' @export liketrace

liketrace = function(x, temp = NULL){

  # x  : data (1 column: temperatures, 2 column likelihood values)
  # temp: temperature to be plotted

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));
  j     = dim(index)[1];

  tempvalue  = round(unique(x$invTemp), 4);

  if(!is.null(temp)){

    ts.plot(x$logL[index[temp, 1]:index[temp, 2]],
            main = paste( "Inverse Temperature ", tempvalue[temp], " - ", temp),
            ylab = "log-likelihood", xlab = "Iteration");

  }else{

    for(i in 1:j){

      aux = index[i,1]:index[i,2];
      ts.plot(x$logL[aux],
              main = paste( "Inverse Temperature ", tempvalue[i], " - ", i),
              ylab = "log-likelihood", xlab = "Iteration");
    }

    cat("Number of traceplots", j);

  }
}
