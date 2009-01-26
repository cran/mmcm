\name{mmcm}
\alias{mmcm}
\docType{package}
\title{
Modified Maximum Contrast Method
}
\description{
This package provides an implementation of the modified maximum contrast method.

Current version supports function \code{\link{mmcm.resamp}} that
gives \eqn{P}-value for the modified maximum contrast statistics
by using a resampling based procedure.

In a one-way problem testing pattern of several factor level means,
maximum contrast statistics may be used.
But under the unbalanced sample size, denominator of the maximum contrast
statistics is overestimated.
Denominator of the modified maximum contrast statistics it is not influenced
by unbalanced sample size.
}
\author{
Author: Kengo NAGASHIMA and Yasunori SATO

Maintainer: Kengo NAGASHIMA \email{nagasima@josai.ac.jp}
}
\references{
Sato, Y., Laird, N.M., Nagashima, K., et al. (2009)
A new statistical screening approach for finding pharmacokinetics-related genes in genome-wide studies.
\emph{Pharmacogenomics J.}, in press.
\url{http://www.ncbi.nlm.nih.gov/pubmed/19104505}
}
\keyword{htest}
%\note{}
\seealso{
  \code{\link{mmcm.resamp}}
}
%\details{
%
%}
