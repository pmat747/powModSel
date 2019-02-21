#' @title Explore a dataset
#' @name expldata
#' @description This function counts the number of samples per temperature.
#' @usage expldata(x)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.
#' @return It produces a data.frame with the number of samples per temperature.
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' expldata(ligoVirgoSim)
#' }
#' @export
expldata = function(x){

  count = stats::aggregate(logL~invTemp, data = x, FUN = length);

  colnames(count) = c("invTemp", "n");

  if(all(diff(count$n) == 0)){

    cat("All temperatures have the same number of samples", "\n")

  }else{

    cat("Minimum number of observations:", "\n")
    print(count[min(count$n) == count$n,])

    cat("Maximum number of observations:", "\n")
    print(count[max(count$n) == count$n,])

  }

  return(count)
}
