\name{print.mmcm.resamp}
\alias{print.mmcm.resamp}
\title{Print function for mmcm.resamp object}
\description{
  This function print result of function \code{\link{mmcm.resamp}}
  (resampling based \eqn{P}-value, maximum contrast and related statistics)
}
\usage{
  \method{print}{mmcm.resamp}(x, ...)
}
\arguments{
  \item{x}{
    Object of class \code{mmcm.resamp}, which is result of
    function \code{\link{mmcm.resamp}}.
  }
  \item{...}{Other parameters used in the function.}
}
\value{
  \item{observed statistics}{
    observed value of the modified maximum contrast statistics (\eqn{T^{\prime}_{\max}}{T'_max}).
  }
  \item{contrast}{
    suffix of coefficient vector and coefficient vector of \eqn{k}th pattern 
    that gives the maximum contrast statistics (colmun number of coefficient matrix).
  }
  \item{P-value}{
    \eqn{P}-value of the modified maximum contrast statistics, which is
    simulated by using a resampling based procedure.
  }

  The case where printed "more than 2 contrast vectors was selected ...",
  some contrast may be unsuitable.
}
\author{
Author: Kengo NAGASHIMA and Yasunori SATO

Maintainer: Kengo NAGASHIMA \email{nagasima@josai.ac.jp}
}
\seealso{
  \code{\link{mmcm.resamp}}
}
\examples{
## Example 1 ##
#  true response pattern: dominant c_3=(1, 1, -2)
set.seed(326584)
x <- c(
  rlnorm(130, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 90, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 10, meanlog = 0.85, sdlog = 0.25)
)
g <- rep(1:3, c(130, 90, 10))

contrast <- cbind(
  c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
)
y <- mmcm.resamp(x, g, contrast, 10000)
print(y)           # print function
y                  #   = method for "mmcm.resamp" class
}
\keyword{print}

