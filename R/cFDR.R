#' Conditional false discovery rate
#' 
#' Estimate the conditional false discovery rate for a set of p-values, 
#' conditional on a second set of p-values
#' 
#' @param data data frame or (numeric) matrix with p-values as columns
#' @param p1 p-values for which the cFDR will be estimated; either the name of the columns in 
#'           \code{data}, or a vectors of p-values
#' @param p2 the p-values on which the cFDR will be conditioned on; either name or vector of p-values
#'           \code{data}, or a vectors of p-values
#' @param p2_threshold cutoff for pre-filtering on \code{p2}: only variants 
#'                     with \code{p2} \leq \code{p2_threshold} are included
#'                     
#' @return A data frame consisting of the p-values and the cFDR
#' 
#' @export
cFDR = function(data, p1, p2, p2_threshold = 1E-3)
{
  "Hello World!"
}
