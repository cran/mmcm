\name{mcm.mvt}
\alias{mcm.mvt}
\title{
  The maximum contrast method by using the quasi-Monte Carlo method
}
\description{
  This function gives \eqn{P}-value for the maximum contrast statistics by using the quasi-Monte Carlo method from \code{\link[mvtnorm:pmvt]{pmvt}} function of package \code{mvtnorm}.
}
\usage{
mcm.mvt(
  x, g, contrast,
  algorithm = GenzBretz()
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values}
  \item{g}{a integer vector giving the group for the corresponding elements of x}
  \item{contrast}{a numeric coefficient matrix for the maximum contrast statistics}
  \item{algorithm}{an object of class '\code{\link[mvtnorm:GenzBretz]{GenzBretz}}' defining the hyper parameters of this algorithm}
}
\details{
  \code{\link{mcm.mvt}} performs the maximum contrast method that is detecting a true response pattern.

  \eqn{Y_{ij} (i=1,~ 2,~\ldots ;~ j=1,~ 2,~ \ldots,~ n_i)}{Y_ij (i = 1, 2, ...; j = 1, 2, ..., n_i)} is an observed response for \eqn{j}th individual in \eqn{i}th group.

  \eqn{\bm{C}}{C} is coefficient matrix for the maximum contrast statistics (\eqn{i \times k}{i x k} matrix, \eqn{i}: No. of groups, \eqn{k}: No. of pattern).
  \deqn{
    \bm{C}=(\bm{c}_1~ \bm{c}_2~ \ldots~ \bm{c}_k)^t
  }{
    C = (c_1 c_2 ... c_k)^t
  }
  \eqn{\bm{c}_k}{c_k} is coefficient vector of \eqn{k}th pattern.
  \deqn{
    \bm{c}_k=(c_{k1}~ c_{k2}~ \ldots~ c_{ki})^t \qquad (\textstyle \sum_i c_{ki}=0)
  }{
    c_k = (c_k1 c_k2 ... c_ki)^t (sum from i of c_ki = 0)
  }
  
  \eqn{T_{\max}}{T_max} is the maximum contrast statistic.
  \deqn{
    \bar{Y}_i=\frac{\sum_{j=1}^{n_i} Y_{ij}}{n_{i}}, ~~
    \bar{\bm{Y}}=(\bar{Y}_1~ \bar{Y}_2~ \ldots~ \bar{Y}_i \ldots~ \bar{Y}_a)^t,
  }{
    Ybar_i = (sum from j of Y_ij) / n_i,
    Ybar = (Ybar_1 Ybar_2 ... Ybar_i ... Ybar_a)^t (a x 1 vector),
  }
  \deqn{
    \bm{D}=diag(n_1,\, n_2,\, \ldots,\, n_i,\, \ldots,\, n_a),~~
    S^2=\frac{1}{\gamma}\sum_{j=1}^{n_i}\sum_{i=1}^{a} (Y_{ij}-\bar{Y}_i)^2,
  }{
    D = diag(n_1, n_2, ..., n_i, ..., n_a) (a x a matrix),
    S^2 = 1/gamma * sum_{j=1}^{n_i} sum_{i=1}^{a} (Y_ij-Ybar_i)^2,
  }
  \deqn{
    \gamma=\sum_{i=1}^{a} (n_i-1),~~
    T_{k}=\frac{\bm{c}^t_k \bar{\bm{Y}}}{\sqrt{S^2\bm{c}^t_k \bm{D} \bm{c}_k}},
  }{
    gamma = sum_{i=1}^{a} (n_i-1),
    T_k = c_k^t Ybar / (S^2 c_k^t D c_k)^(1/2),
  }
  \deqn{
    T_{\max}=\max(T_1,~ T_2,~ \ldots,~ T_k).
  }{
    T_max = max(T_1, T_2, ..., T_k).
  }
  
  Consider testing the overall null hypothesis \eqn{H_0: \mu_1=\mu_2=\ldots=\mu_i}, versus alternative hypotheses \eqn{H_1} for response petterns  (\eqn{H_1: \mu_1<\mu_2<\ldots<\mu_i,~ \mu_1=\mu_2<\ldots<\mu_i,~ \mu_1<\mu_2<\ldots=\mu_i}).
  The \eqn{P}-value for the probability distribution of \eqn{T_{\max}}{T_max} under the overall null hypothesis is
  \deqn{
    P\mbox{-value}=\Pr(T_{\max}>t^{*}_{\max} \mid H_0)
  }{
    P-value = Pr(T_max > t*_max | H0)
  }
  \eqn{t^{*}_{\max}}{t*_max} is observed value of statistics.
  This function gives distribution of \eqn{T_{\max}}{T_max} by using the quasi-Monte Carlo method from package \code{mvtnorm}.
}
\value{
  \item{cont}{
    suffix of coefficient vector of \eqn{k}th pattern that gives maximum contrast statistics (colmun number of coefficient matrix).
  }
  \item{pval}{
    \eqn{P}-value of the maximum contrast statistics, which is simulated by using the quasi-Monte Carlo method.
  }
  \item{astat}{
    vector of observed values (\eqn{T_k}{T_k}) for each the contrast statistics.
  }
  \item{acont}{
    coefficient matrix (copy of argument).
  }
  return value is \code{mmcm.res} class object.
}
%\note{
% x list version
% g factor object version
%}
\references{
Yoshimura, I., Wakana, A., Hamada, C. (1997).
A performance comparison of maximum contrast methods to detect dose dependency.
\emph{Drug Information J.} \bold{31}: 423--432.
}
\author{
Author: Kengo NAGASHIMA and Yasunori SATO

Maintainer: Kengo NAGASHIMA \email{nagasima@josai.ac.jp}
}
\seealso{
  \code{\link[mvtnorm:pmvt]{pmvt}},
  \code{\link[mvtnorm:GenzBretz]{GenzBretz}},
  \code{\link{mmcm.mvt}}
}
\examples{
## Example 1 ##
#  true response pattern: dominant model c=(1, 1, -2)
set.seed(326584)
x <- c(
  rlnorm(130, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 90, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 10, meanlog = 0.85, sdlog = 0.25)
)
g <- rep(1:3, c(130, 90, 10))
boxplot(
  x ~ g,
  width=c(length(g[g==1]),length(g[g==2]),
    length(g[g==3])),
  main="Dominant model (sample data)",
  xlab="Genotype", ylab="PK parameter"
)

# coefficient matrix
# c_1: additive, c_2: recessive, c_3: dominant
contrast <- rbind(
  c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
)
y <- mcm.mvt(x, g, contrast)
y

## Example 2 ##
#  for dataframe
#  true response pattern: pos = 1 dominant  model c=( 1,  1, -2)
#                               2 additive  model c=(-1,  0,  1)
#                               3 recessive model c=( 2, -1, -1)
set.seed(8415849)
x   <- c(
  rlnorm(130, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 90, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 10, meanlog = 0.85, sdlog = 0.25),
  rlnorm(130, meanlog = 0.79, sdlog = 0.1),
  rlnorm( 90, meanlog = 0.85, sdlog = 0.1),
  rlnorm( 10, meanlog = 0.91, sdlog = 0.25),
  rlnorm(130, meanlog = 0.85, sdlog = 0.1),
  rlnorm( 90, meanlog = 0.91, sdlog = 0.1),
  rlnorm( 10, meanlog = 0.91, sdlog = 0.25)
)
g   <- rep(rep(1:3, c(130, 90, 10)), 3)
pos <- rep(c("rsXXXX", "rsYYYY", "rsZZZZ"), each=230)
xx <- data.frame(pos = pos, x = x, g = g)

# coefficient matrix
# c_1: additive, c_2: recessive, c_3: dominant
contrast <- rbind(
  c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
)
mmcmtapply <- function(r) {
  mcm.mvt(
    xx$x[xx$pos==r[1]], xx$g[xx$pos==r[1]],
    contrast
  )
}
y <- tapply(xx$pos, xx$pos, mmcmtapply)
yy <- data.frame(
  Pos = as.vector(names(y)),
  Pval = as.vector(sapply(y, "[[", 3)),
  Pattern = as.vector(sapply(y, "[[", 4))
)
# miss-detection!
yy
}
\keyword{htest}
