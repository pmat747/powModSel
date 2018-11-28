#' @title Log-likelihood trace plot
#' @description This function produces the traceplot for the log-likelihood at each temperature
#' @param x dataframe
#' @return A set of traceplots
#' @export liketrace

liketrace = function(x, num = NULL){

  # x  : data (1 column: temperatures, 2 column likelihood values)
  # num: temperature to be plotted

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));
  j     = dim(index)[1];

  temp  = round(unique(x$invTemp), 4);

  if(!is.null(num)){

    ts.plot(x$logL[index[num, 1]:index[num, 2]],
            main = paste( "Inverse Temperature ", temp[num], " - ", num),
            ylab = "log-likelihood", xlab = "Iteration");

  }else{

    for(i in 1:j){

      aux = index[i,1]:index[i,2];
      ts.plot(x$logL[aux],
              main = paste( "Inverse Temperature ", temp[i], " - ", i),
              ylab = "log-likelihood", xlab = "Iteration");
    }

    cat("Number of traceplots", j);

  }
}
