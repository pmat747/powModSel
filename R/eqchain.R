#' @title equal chain lengths
#' @description This function equates the Markov chain lengths at all temperatures according to the shortest one
#' @param x dataframe
#' @return A numeric value, and an optional plot
#' @export eqchain

eqchain = function(x){

  count = stats::aggregate(logL~invTemp, data = x, FUN = length);
  count = count$logL; # Number of elements per temperature

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
  ### ###
  return(x);
}
