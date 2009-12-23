\name{print.mmcm.res}
\alias{print.mmcm.res}
\title{Print function for mmcm.res object}
\description{
  This function print result of function \code{\link{mmcm.mvt}} and
  \code{\link{mmcm.resamp}}
}
\usage{
  \method{print}{mmcm.res}(x, ...)
}
\arguments{
  \item{x}{
    Object of class \code{mmcm.res}, which is result of
    function \code{\link{mmcm.mvt}}, \code{\link{mmcm.resamp}} and \code{\link{mcm.mvt}}.
  }
  \item{...}{Other parameters used in the function.}
}
\value{
  \item{observed statistics}{
    observed value of the maximum contrast statistics / a modified maximum contrast statistics (\eqn{T_{\max}}{T_max} / \eqn{T^{\prime}_{\max}}{T'_max} / \eqn{T^{\prime\prime}_{\max}}{T''_max}).
  }
  \item{contrast}{
    suffix of coefficient vector and coefficient vector of \eqn{k}th pattern that gives the maximum contrast statistics / a modified maximum contrast statistics (colmun number of coefficient matrix).
  }
  \item{P-value}{
    \eqn{P}-value of the maximum contrast statistics / a modified maximum contrast statistics
  }

  The case where printed "more than 2 contrast vectors was selected ...", some contrast may be unsuitable.
}
\author{
Author: Kengo NAGASHIMA and Yasunori SATO

Maintainer: Kengo NAGASHIMA \email{nagasima@josai.ac.jp}
}
\seealso{
  \code{\link{mmcm.mvt}}, \code{\link{mmcm.resamp}}, \code{\link{mcm.mvt}}
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

contrast <- rbind(
  c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
)
y <- mmcm.resamp(x, g, contrast, 10000)
print(y)           # print function
y                  #   = S3 method for "mmcm.res" class
y <- mmcm.mvt(x, g, contrast)
print(y)           # print function
y                  #   = S3 method for "mmcm.res" class
y <- mcm.mvt(x, g, contrast)
print(y)           # print function
y                  #   = S3 method for "mmcm.res" class
}
\keyword{print}

