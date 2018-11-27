#' @title Auxiliary functions
#' @description These functions are transversally used by the rest of functions
#' @param
#' @return
#' @export

#######################################
### log(x+y)=logplus(log(x),log(y)) ###
#######################################

logplus <- function(x,y)
{
  if(x>y) x + log(1+exp(y-x))
  else    y + log(1+exp(x-y))
}

# logplus function for a vector x
logplusvec = function(x){
  r = -Inf;
  for(i in x){
    r = logplus(r, i);
  }
  return(r);
}
