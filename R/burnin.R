#' @title Burnin period
#' @description It deletes initial likelihood values at each temperature as burnin period
#' @param x dataframe
#' @param burnin integer which stands for the number of values to be deleted
#' @param percentage percentage of the burnin period
#' @return A numeric value, and an optional plot
#' @return A dataframe
#' @export burnin

burnin = function(x, n_burnin, percentage = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));

  if(!is.null(percentage)){

    if((percentage < 0) && (percentage > 100)){
      stop("Percentage must be in [0-100]")}

    n_burnin= round((index[,2] - index[,1] + 1) * percentage / 100);

  }else if((n_burnin %% 1 != 0) || (n_burnin < 0)){
    stop("n_burnin must be a non-negative integer")
  }

  index[, 1] = index[, 1] + n_burnin;

  index = unlist(apply(index, 1, function(y)seq(y[1], y[2], by=1)) );

  #out   = x[as.character(index),];
  out   = x[index,];
  rownames(out) = NULL; # reseting rownames # need when taking row values
  out   = as.data.frame(out);

  return(out);

}
