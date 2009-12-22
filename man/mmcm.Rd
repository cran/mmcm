\name{mmcm}
\alias{mmcm}
\docType{package}
\title{
Modified Maximum Contrast Method
}
\description{
This package provides an implementation of a modified maximum contrast method and the maximum contrast method.

This version supports functions \code{\link{mmcm.mvt}}, \code{\link{mcm.mvt}} that gives P-value by using the quasi-Monte Carlo method, and \code{\link{mmcm.resamp}} that gives P-value by using the permutation method.

In a one-way problem testing pattern of several factor level means, the maximum contrast statistics (Yoshimura, I., 1997) may be used. But under the unequal sample size situation, denominator of the maximum contrast statistics is overestimated. Thus we propose a modified maximum contrast statistics for the unequal sample size situation. Denominetor of a modified maximum contrast statistics is not influenced under the unequal sample size situation.
}
\author{
Author: Kengo NAGASHIMA and Yasunori SATO

Maintainer: Kengo NAGASHIMA \email{nagasima@josai.ac.jp}
}
\references{
Sato, Y., Laird, N.M., Nagashima, K., et al. (2009).
A new statistical screening approach for finding pharmacokinetics-related genes in genome-wide studies.
\emph{Pharmacogenomics J.} \bold{9}(2): 137--146.
\url{http://www.ncbi.nlm.nih.gov/pubmed/19104505}

Yoshimura, I., Wakana, A., Hamada, C. (1997).
A performance comparison of maximum contrast methods to detect dose dependency.
\emph{Drug Information J.} \bold{31}: 423--432.
}
\keyword{htest}
%\note{}
\seealso{
  \code{\link{mmcm.mvt}}, \code{\link{mmcm.resamp}}, \code{\link{mcm.mvt}}
}
%\details{
%
%}
