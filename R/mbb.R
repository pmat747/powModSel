#' @title Moving block bootstrap
#' @description
#' @param data A dataframe
#' @return list
#' @export mbb

mbb = function(x, l_mbb, n_bootstrap, temps = NULL){

  count = aggregate(logL~invTemp, data = x, FUN = length);# (temp, length)
  count = count$logL; # Number of elements per temperature
  # unique(x$temperature); # temperatures in x
  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));# indexes

  if( !is.null(temps)){ # selecting certain temperatures (temps)
    index = index[temps, ]; # for original dataset
    count = count[temps];
    newX  = unlist(apply(index, 1, function(x)seq(x[1],x[2], by = 1)));
    newX  = x[as.character(newX),];
    rownames(newX) = NULL; # reseting rownames
    x      = newX;  # Copied to be used below
    ref_ti = ti(newX); # true values for the subset
    ref_ss = ss(newX);

    index = which(diff(x$invTemp)!=0); # for new dataset
    index = cbind(c(1, index + 1), c(index, dim(x)[1]));# indexes

  }else{

    ref_ti = ti(x); # estimate using original data
    ref_ss = ss(x);

  }

  b_mbb = ceiling(count / l_mbb);# Numbers of blocks to be sampled per temperature (rounded up)
  minb_mbb = min(b_mbb); # minimum (b_mbb * l_mbb >= count)
  N_mbb = count - l_mbb + 1; # Number of blocks per temperature
  minCount = min(count); # Minimum chain length.
  minN_mbb = min(N_mbb); # Sampling according to the minimum number of blocks
  n_temp = length(count); # Number of temperatures
  Zs      = NULL;
  indx   = c(apply(matrix(index[,1] - 1), 1, function(x) rep(x, minCount)));

  Sam  = sample(x = minN_mbb, size = n_bootstrap * minb_mbb, replace = TRUE);# Sampling all the Temps at the same position
  # Each row contains starting block positions to generate a bootstrap dataset
  Sam  = matrix(Sam, nrow = n_bootstrap);

  ptime = proc.time()[1];

  for(bi in 1:n_bootstrap){

    if (bi %% 10 == 0){
      print(paste("Iteration", bi, ",", "Time elapsed",
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes"));
    }

    # Boostrap data
    sam = bsam(Sam[bi, ], l_mbb, indx, minCount, n_temp);
    x1  = x[as.character(sam), ];
    rownames(x1) = NULL; # reseting rownames

    Zs  = rbind(Zs, c(ti(x1), ss(x1)));

  } # end bootstrap

  colnames(Zs) = c("ti","ss")

  res = cbind(Zs[,1] - ref_ti, Zs[,2] - ref_ss);
  colnames(res) = c("ti","ss");

  return(list(Zs = Zs, se = apply(Zs, 2, sd), res = res) );

}

bsam = function(y, l_mbb, indx, minCount, n_temp){

  # Internally used in 'mbb' function

  # This function takes the starting block positions.
  # It produces the block positions for each temperature.
  # Inputs:
  # - y: starting block positions (first temperature);
  # - l_mbb: block length;
  # - indx: to calculate positions in other temperatures;
  # - minCount: minimum chain length;
  # - n_temp: number of temperatures;

  out = c(apply(matrix(y), 1, function(x) seq(from = x, length = l_mbb)));
  out = out[1:minCount]; # discarding overleft points
  out = rep(out, n_temp); # No need because recycle rule
  #indx = c(apply(matrix(index[,1] - 1), 1, function(x) rep(x, minCount))); # Outside
  return(out + indx);

}

