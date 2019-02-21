#' @title Equal chain lengths
#' @name eqchain
#' @description This function equates the Markov chain lengths at all the temperatures according to the shortest one or a specific value.
#' @usage eqchain(x,n)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.
#' @param n Sequence length per temperature.
#' @return It produces a new \code{x} with a equal number of log-likelihood values at each temperature.
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' eqchain(ligoVirgoSim)
#' }
#' @export
eqchain = function(x, n = NULL){

  count = stats::aggregate(logL~invTemp, data = x, FUN = length);
  count = count$logL; # Number of elements per temperature

  if( is.null(n) ){

    if(all(diff(count) == 0)){
      cat("All Markov chains have same lengths", "\n")
      return(x)
    }

    minCount = min(count);

    index = which(diff(x$invTemp)!=0);
    index = cbind(c(1, index + 1), c(index, dim(x)[1]));

    aux   = apply(matrix(index[,1]), 1,
                  function(x) seq(from = x, length = minCount));

    aux   = c(aux);

    x     = x[aux, ];
    rownames(x) = NULL; # reseting rownames
    x     = as.data.frame(x);

    ### Print on screen ###
    r = as.data.frame(cbind(1:length(count), abs(count - minCount)));
    colnames(r) = c("invTemp", "Deleted");
    cat("Number of observations (from the tail) deleted per temperature", "\n");
    print(r);
    cat("New chain lengths:", minCount,"\n");

    return(x);

  }else{

    minCount = min(count);
    if(n > minCount) stop(paste("n must be lower than ", minCount,
                                ", the smallest group", sep = ""));
    # Same as above
    index = which(diff(x$invTemp)!=0);
    index = cbind(c(1, index + 1), c(index, dim(x)[1]));
    aux   = apply(matrix(index[,1]), 1,
                  function(x) seq(from = x, length = n));
    aux   = c(aux);
    x     = x[aux, ];
    rownames(x) = NULL; # reseting rownames
    x     = as.data.frame(x);
    return(x);
  }
}
