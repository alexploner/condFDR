% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/condFDR.R
\name{cFDR}
\alias{cFDR}
\title{Conditional false discovery rate}
\usage{
cFDR(data, p1, p2, p2_threshold = 0.001, mc.cores = 1)
}
\arguments{
\item{data}{data frame or (numeric) matrix with p-values as columns}

\item{p1}{p-values for which the cFDR will be estimated; either the name/index of the columns in
\code{data}, or a vectors of p-values}

\item{p2}{the p-values on which the cFDR will be conditioned on; either name or vector of p-values
\code{data}, or a vectors of p-values}

\item{p2_threshold}{cutoff for pre-filtering on \code{p2}: only variants
with \code{p2} \eqn{\leq}{<=} \code{p2_threshold} are included}

\item{mc.cores}{number of cores to use for parallel calculation; defaults to one
(i.e. no parallel calculation), but should absolutely be increased if
your system supports it, as this will speed up execution very nicely.}
}
\value{
A data frame with column cFDR: if \code{data} was specified, the column
is simply added at the end; if only \code{p1} and \code{p2} were specified, a
data frame with three columns (the original p-values and the cFDR).
}
\description{
Estimate the conditional false discovery rate for a set of p-values,
conditional on a second set of p-values
}
\examples{
data(psynth)
res = cFDR(psynth, "p1", "p2", p2_threshold = 1E-5)
head(res)
head(subset(res, cFDR < 0.01))
}
\seealso{
\code{\link[parallel]{mclapply}} for details on parallel calculations and \code{mc.cores}
}
