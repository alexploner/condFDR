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
#'                     with \code{p2} \leq \code{p2_threshold} are included
#'
#' @return A data frame with column cFDR: if \code{data} was specified, the column
#' is simply added at the end; if only \code{p1} and \code{p2} were specified, a
#' data frame with three columns (the original p-values and the cFDR).
#'
#' @export
cFDR = function(data, p1, p2, p2_threshold = 1E-3)
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

  ## Loop
  nn    <- length(p1)
  denom <- rep(0, nn)
  for (i in 1:nn) {

    ## The edge point
    x <- p1[i]
    y <- p2[i]

    ## Vectors
    dd = p2 <= y
    ee = dd & (p1 <= x)

    ## Combine
    denom[i] = length(which(ee)) / length(which(dd))

  }

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
