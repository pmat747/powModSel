#' @title Moving Block Bootstrap
#' @name mbb
#' @description This function produces a set of marginal likelihood estimates for moving block bootstrap observations.  Its main use is for calculating the standard error associated to the thermodynamic integration and stepping-stone sampling estimates.
#' @usage mbb(x, bl, nboot, temp = NULL)
#' @param x A data frame with the folloging columns: \code{logL} and \code{invTemp}, which contain the log-likelihood and inverse temperature values, respectively.
#' @param bl Block lenghts.
#' @param nboot Number of bootstrap observations to be analysed.
#' @param temp It indicates the temperatures to be used in the analysis, for instance, c(1,3,K) considers the temperatures at those positions, where \eqn{K} is the number of temperatures.  In this case, the temperatures must be sorted in an increasing order.  Note that samples from the prior and posterior must be included in the process.  \code{NULL} stands for all the temperatures in \code{x}.
#' @details For a block length equal to 1 (\code{bl}=1) the original bootstrap method for i.i.d. data is recovered.  A block lenght greater than one allows to take into account a potential autocorrelation within the Markov chains.  \code{mbb} is being designed to take also into account potential cross-correlation between the Markov chains due to swaps in parallel tempering sampling.  See more details in Maturana R. et al. (2018)
#' @return A list containing the following components:
#'         \item{Zs}{Marginal likelihood estimates via \code{ti} and \code{ss}.}
#'         \item{se}{Standard deviation of the marginal likelihood estimates calculated for the bootstrap observations.}
#'         \item{res}{Marginal likelihood estimate differences between the ones calculated for the bootstrap observations and the original dataset \code{x}.}
#' @author Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
#' @references Kunsch, H. R. 1989. The Jackknife and the Bootstrap for General Stationary Observations. \emph{The Annals of Statistics} \bold{17}(3), 1217--1241.
#' Maturana Russel, P., Meyer, R., Veitch, J., and Christensen, N. 2018. The stepping-stone algorithm for calculating the evidence of gravitational wave models. arXiv preprint arXiv:1810.04488
#' @examples
#' \dontrun{
#' data(ligoVirgoSim)
#' R = mbb(ligoVirgoSim, bl = 10, nboot = 20, temp = NULL)
#' R$se; # standard error of the marginal likelihood estimates
#' }
#' @export
mbb = function(x, bl, nboot, temp = NULL){

  count = stats::aggregate(logL~invTemp, data = x, FUN = length);# (temp, length)
  count = count$logL; # Number of elements per temperature
  # unique(x$temperature); # temperatures in x
  index = which(diff(x$invTemp)!=0);
  index = cbind(c(1, index + 1), c(index, dim(x)[1]));# indexes

  if( !is.null(temp)){ # selecting certain temperatures (temp)
    index = index[temp, ]; # for original dataset
    count = count[temp];
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

  b_mbb = ceiling(count / bl);# Numbers of blocks to be sampled per temperature (rounded up)
  minb_mbb = min(b_mbb); # minimum (b_mbb * bl >= count)
  N_mbb = count - bl + 1; # Number of blocks per temperature
  minCount = min(count); # Minimum chain length.
  minN_mbb = min(N_mbb); # Sampling according to the minimum number of blocks
  n_temp = length(count); # Number of temperatures
  Zs      = NULL;
  indx   = c(apply(matrix(index[,1] - 1), 1, function(x) rep(x, minCount)));

  Sam  = sample(x = minN_mbb, size = nboot * minb_mbb, replace = TRUE);# Sampling all the invTemp at the same position
  # Each row contains starting block positions to generate a bootstrap dataset
  Sam  = matrix(Sam, nrow = nboot);

  ptime = proc.time()[1];

  for(bi in 1:nboot){

    if (bi %% 10 == 0){
      print(paste("Iteration", bi, ",", "Time elapsed",
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes"));
    }

    # Boostrap data
    sam = bsam(Sam[bi, ], bl, indx, minCount, n_temp);
    x1  = x[as.character(sam), ];
    rownames(x1) = NULL; # reseting rownames

    Zs  = rbind(Zs, c(ti(x1), ss(x1)));

  } # end bootstrap

  colnames(Zs) = c("ti","ss")

  res = cbind(Zs[,1] - ref_ti, Zs[,2] - ref_ss);
  colnames(res) = c("ti","ss");

  return(list(Zs = Zs, se = apply(Zs, 2, stats::sd), res = res) );

}

bsam = function(y, bl, indx, minCount, n_temp){

  # Internally used in 'mbb' function

  # This function takes the starting block positions.
  # It produces the block positions for each temperature.
  # Inputs:
  # - y: starting block positions (first temperature);
  # - bl: block length;
  # - indx: to calculate positions in other temperatures;
  # - minCount: minimum chain length;
  # - n_temp: number of temperatures;

  out = c(apply(matrix(y), 1, function(x) seq(from = x, length = bl)));
  out = out[1:minCount]; # discarding overleft points
  out = rep(out, n_temp); # No need because recycle rule
  #indx = c(apply(matrix(index[,1] - 1), 1, function(x) rep(x, minCount))); # Outside
  return(out + indx);

}

