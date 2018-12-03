#' @title Burnin period
#' @description It deletes initial likelihood values at each temperature as burnin period
#' @param x dataframe
#' @param percentage percentage of the burnin period
#' @param burnin integer which stands for the number of values to be deleted
#' @return A dataframe
#' @export burnin

burnin = function(x, percentage = NULL, nburnin = NULL){

  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));

  if(!is.null(percentage)){

    if((percentage < 0) || (percentage > 100)){
      stop("Percentage must be in [0-100]")}

    nburnin= round((index[,2] - index[,1] + 1) * percentage / 100);

  }else if(!is.null(nburnin)){

    if((nburnin %% 1 != 0) || (nburnin < 0)){
      stop("nburnin must be a non-negative integer")}

    aux = min(index[,2] - index[,1]);

    if( aux < nburnin){
      stop(paste("nburnin must be lowest than ", aux, sep = "", "\n"))
    }

  }else{
    stop("Burnin period must be specified")
  }

  index[, 1] = index[, 1] + nburnin;

  index = unlist(apply(index, 1, function(y)seq(y[1], y[2], by=1)) );

  #out   = x[as.character(index),];
  out   = x[index,];
  rownames(out) = NULL; # reseting rownames # need when taking row values
  out   = as.data.frame(out);

  return(out);

}
