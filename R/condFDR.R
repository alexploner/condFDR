#' Conditional false discovery rate
#'
#' Estimate the conditional false discovery rate for a set of p-values,
#' conditional on a second set of p-values
#'
#' @param data data frame or (numeric) matrix with p-values as columns
#' @param p1 p-values for which the cFDR will be estimated; either the name/index of the columns in
#'           \code{data}, or a vectors of p-values
#' @param p2 the p-values on which the cFDR will be conditioned on; either name or vector of p-values
#'           \code{data}, or a vectors of p-values
#' @param p2_threshold cutoff for pre-filtering on \code{p2}: only variants
#'                     with \code{p2} \eqn{\leq}{<=} \code{p2_threshold} are included
#' @param mc.cores number of cores to use for parallel calculation; defaults to one
#'                 (i.e. no parallel calculation), but should absolutely be increased if
#'                 your system supports it, as this will speed up execution very nicely.
#'
#' @return A data frame with column cFDR: if \code{data} was specified, the column
#' is simply added at the end; if only \code{p1} and \code{p2} were specified, a
#' data frame with three columns (the original p-values and the cFDR).
#'
#' @export
#' @seealso \code{\link[parallel]{mclapply}} for details on parallel calculations and \code{mc.cores}
#' @examples
#' data(psynth)
#' res = cFDR(psynth, "p1", "p2", p2_threshold = 1E-5)
#' head(res)
#' head(subset(res, cFDR < 0.01))
cFDR = function(data, p1, p2, p2_threshold = 1E-3, mc.cores = 1)
{
  ## Extract the data
  if ( !missing(data) ) {
    stopifnot( is.matrix(data) | is.data.frame(data)  )
    p1 = data[, p1]
    p2 = data[, p2]
  } else {
    stopifnot( length(p1) == length(p2))
    p1_name = deparse(substitute(p1))
    p2_name = deparse(substitute(p2))
  }

  ## Check: probabilities, no missing values
  doCheck = function(x) is.numeric(x) & !any(is.na(x)) & !any(x<0) & !any(x>1)
  stopifnot( doCheck(p1) & doCheck(p2) )

  ## Subset
  ndx = p2 <= p2_threshold
  stopifnot( any(ndx) )
  p1  = p1[ndx]
  p2  = p2[ndx]

  ## Loop: unwrapped in mcapply
  nn = length(p1)
  calc.denom = function(i)
  {
    ## The edge point
    x = p1[i]
    y = p2[i]

    ## Vectors
    dd = p2 <= y
    ee = dd & (p1 <= x)

    ## Combine
    sum(ee) / sum(dd)
  }
  denom = unlist( parallel::mclapply(1:nn, calc.denom, mc.cores = mc.cores) )

  ## Calibrate the p-values
  cfdr = p1 / denom

  ## Build output
  if (!missing(data) ) {
    ret = cbind(data[ndx, ], cFDR = cfdr)
  } else {
    ret = data.frame(p1, p2, cfdr)
    colnames(ret) = c(p1_name, p2_name, "cFDR")
  }

  ret
}


#' Conjunctional conditional false discovery rate
#'
#' Estimate the conjunctional conditional false discovery rate for two set of p-values,
#' conditional on each other, as measure for pleiotropy.
#'
#' @param data data frame or (numeric) matrix with p-values as columns
#' @param p1 p-values for which the ccFDR will be estimated; either the name/index of the columns in
#'           \code{data}, or a vectors of p-values
#' @param p2 the p-values on which the ccFDR will be conditioned on; either name or vector of p-values
#'           \code{data}, or a vectors of p-values
#' @param p_threshold cutoff for pre-filtering the p-values; a vector of either length one (in which
#'          case both sets of pvalues have the same threshold) or length two (where the thresholds refer
#'          to \code{p1} and \code{p2} in that order).
#' @param mc.cores number of cores to use for parallel calculation; defaults to one
#'                 (i.e. no parallel calculation), but should absolutely be increased if
#'                 your system supports it, as this will speed up execution very nicely.
#'
#' @return A data frame with columns cFDR1, cFDR2, and ccFDR: if \code{data} was specified,
#' the columns are simply added at the end; if only \code{p1} and \code{p2} were specified, a
#' data frame with five columns (the original p-values, plus the cFDRs and the ccFDR).
#'
#' @export
#' @seealso \code{\link[parallel]{mclapply}} for details on parallel calculations and \code{mc.cores}
#' @examples
#' data(psynth)
#' res = ccFDR(psynth, "p1", "p2", p_threshold = 1E-5)
#' head(res)
#' head(subset(res, ccFDR < 0.01))
ccFDR = function(data, p1, p2, p_threshold = 1E-3, mc.cores = 1)
{
  ## Extract the data
  if ( !missing(data) ) {
    stopifnot( is.matrix(data) | is.data.frame(data)  )
    p1 = data[, p1]
    p2 = data[, p2]
  } else {
    stopifnot( length(p1) == length(p2))
    p1_name = deparse(substitute(p1))
    p2_name = deparse(substitute(p2))
  }

  ## Check: probabilities, no missing values
  doCheck = function(x) is.numeric(x) & !any(is.na(x)) & !any(x<0) & !any(x>1)
  stopifnot( doCheck(p1) & doCheck(p2) )

  ## Subset
  p_threshold = rep_len(p_threshold, length.out = 2)
  ndx = (p1 <= p_threshold[1]) | (p2 <= p_threshold[2])
  stopifnot( any(ndx) )
  p1  = p1[ndx]
  p2  = p2[ndx]

  ## Loop
  nn = length(p1)
  calc.denoms = function(i)
  {
    ## The edge point
    x = p1[i]
    y = p2[i]

    ## Vectors
    dd1 = p2 <= y
    dd2 = p1 <= x
    ee = dd1 & dd2

    ## Combine
    ee_n = length(which(ee))
    denom1 =  ee_n / length(which(dd2))
    denom2 =  ee_n / length(which(dd2))

    c(denom1, denom2)
  }
  ## This needs a bit care: unlist generates a vector with alternating denom1/denom2
  ## We can pour this vector into a matrix with two rows, without having to re-order,
  ## and can extract the rows in the next step
  denoms = matrix( unlist( parallel::mclapply(1:nn, calc.denoms, mc.cores = mc.cores) ), nrow = 2)

  ## Calibrate the p-values; note: use rows!
  cfdr1 = p1 / denoms[1, ]
  cfdr2 = p2 / denoms[2, ]
  ccfdr = pmax(cfdr1, cfdr2)

  ## Build output
  if ( !missing(data) ) {
    ret = cbind(data[ndx, ], cFDR1 = cfdr1, cFDR2 = cfdr2, ccFDR = ccfdr)
  } else {
    ret = data.frame(p1, p2, cFDR1 = cfdr1, cFDR2 = cfdr2, ccFDR = ccfdr)
    colnames(ret)[1:2] = c(p1_name, p2_name)
  }

  ret
}


#' Synthetic SNP p-values
#'
#' Synthetic p-values for two phenotypes across a moderately larger number of
#' (fake) SNPs.
#'
#' @format A data frame with 1,750,000 rows and 3 variables:
#' \describe{
#'   \item{SNP}{ids of fictitious SNPs}
#'   \item{p1}{p-values for first phenotype}
#'   \item{p2}{p-values for second phenotype}
#' }
"psynth"
