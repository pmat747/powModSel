#' @title Stepping-stone sampling
#' @description It calculates a stepping-stone sampling estimate given the likelihood and the temperature values
#' @param x A dataframe
#' @param temp Integer vector subset of 1:K, where K:number of temperatures
#' @return Stepping-stone sampling estimate (numeric value)
#' @export ss

ss = function(x, temp = NULL){

  if( !is.null(temp)){ # selecting certain temperatures (temps)
    #count = stats::aggregate(logL~invTemp, data = x, FUN = length);
    index = which(diff(x$invTemp)!=0);
    index = cbind(c(1, index + 1), c(index, dim(x)[1]));
    index = index[temp,];
    newX  = unlist(apply(index, 1, function(x)seq(x[1],x[2], by = 1)));
    newX  = x[as.character(newX),];
    #rownames(newX) = NULL; # reseting rownames
    x = newX;
  }

  rownames(x) = NULL; # reseting rownames

  count = stats::aggregate(logL~invTemp, data = x, FUN = length);
  # 'count' could be sorted here, increasing order, if 'x' is not ordered
  diff  = diff(count$invTemp); #=diff(unique(x$temperature));#difference between temperatures
  count = count$logL[-dim(count)[1]]; # Number of elements per temperature/no posterior
  delta = rep(diff, times = count); # replicating diff in original data
  N     = sum(count); # = length(delta)
  logL  = x[1:N, 'logL']; # extracting loglike
  loglDelta = logL * delta; # loglike x delta
  invTemp   = x[1:N, 'invTemp']; # extracting temperature
  Rss   = cbind(invTemp, loglDelta);
  Rss   = stats::aggregate(loglDelta~invTemp, data = Rss, FUN = logplusvec);
  Rss   = Rss["loglDelta"] - log(count); # Individual rate estimates

  return(sum(Rss));

}
