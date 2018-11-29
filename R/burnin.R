#' @title Burnin period
#' @description It deletes initial likelihood values at each temperature as burnin period
#' @param x dataframe
#' @param percentage percentage of the burnin period
#' @param burnin integer which stands for the number of values to be deleted
#' @return A dataframe
#' @export burnin

burnin = function(x, percentage = NULL, n_burnin = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));

  if(!is.null(percentage)){

    if((percentage < 0) || (percentage > 100)){
      stop("Percentage must be in [0-100]")}

    n_burnin= round((index[,2] - index[,1] + 1) * percentage / 100);

  }else if(!is.null(n_burnin)){

    if((n_burnin %% 1 != 0) || (n_burnin < 0)){
      stop("n_burnin must be a non-negative integer")}

    aux = min(index[,2] - index[,1]);

    if( aux < n_burnin){
      stop(paste("n_burnin must be lowest than ", aux, sep = "", "\n"))
    }

  }else{
    stop("Burnin period must be specified")
  }

  index[, 1] = index[, 1] + n_burnin;

  index = unlist(apply(index, 1, function(y)seq(y[1], y[2], by=1)) );

  #out   = x[as.character(index),];
  out   = x[index,];
  rownames(out) = NULL; # reseting rownames # need when taking row values
  out   = as.data.frame(out);

  return(out);

}
