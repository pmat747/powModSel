#' @title Burnin period
#' @name burnin
#' @description This function discards the first log-likelihood values at each inverse temperature as burnin period.
#' @usage burnin(x, percentage = NULL, nburnin = NULL)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.
#' @param percentage Percentage of log-likelihod values to discard at each temperature.
#' @param nburnin Number of log-likelihod values to discard at each temperature.
#' @return It produces a new \code{x} with a burnin period.
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' burnin(ligoVirgoSim, percentage=10); #discarding 10\% of log-likelihood values at each temperature
#' }
#' @export
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
