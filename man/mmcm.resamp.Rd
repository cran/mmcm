\name{mmcm.resamp}
\alias{mmcm.resamp}
\title{
  P-value for the modified maximum contrast
  statistics by using a resampling based procedure
}
\description{
  This function gives \eqn{P}-value for the modified maximum
  contrast statistics by using a resampling based procedure.
}
\usage{
mmcm.resamp(
  x, g, contrast,
  nsample = 20000, seed = unclass(Sys.time())
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data values}
  \item{g}{a integer vector giving the group for the corresponding elements of x}
  \item{contrast}{a numeric coefficient matrix for the modified maximum contrast statistics}
  \item{nsample}{specifies the number of resamples (defalt: 20,000)}
  \item{seed}{a 32-bit integer (\eqn{-2^{31}+1\leq \mbox{seed} \leq 2^{31}-1}{-2^31+1 <= seed <= 2^31-1})
    used as the seed for the pseudo-random number generator
    used for resampling (default: time-dependent integer)}
}
\details{
  \code{\link{mmcm.resamp}} performs the modified maximum contrast method that is
  detecting a true response pattern when unbalanced sample size.

  \eqn{Y_{ij} (i=1,~ 2,~\ldots ;~ j=1,~ 2,~ \ldots,~ n_i)}{Y_ij (i = 1, 2, ...; j = 1, 2, ..., n_i)}
  is an observed response for \eqn{j}th individual in \eqn{i}th group.

  \eqn{\bm{C}}{C} is coefficient matrix for the modified maximum contrast statistics 
  (\eqn{i \times k}{i x k} matrix, \eqn{i}: No. of groups, \eqn{k}: No. of pattern).
  \deqn{
    \bm{C}=(\bm{c}_1~ \bm{c}_2~ \ldots~ \bm{c}_k)
  }{
    C = (c_1 c_2 ... c_k)
  }
  \eqn{\bm{c}_k}{c_k} is coefficient vector of \eqn{k}th pattern.
  \deqn{
    \bm{c}_k=(c_{k1}~ c_{k2}~ \ldots~ c_{ki})^t \qquad (\textstyle \sum_i c_{ki}=0)
  }{
    c_k = (c_k1 c_k2 ... c_ki)^t (sum from i of c_ki = 0)
  }
  
  \eqn{T^{\prime}_{\max}}{T'_max} is the modified maximum contrast statistics. 
  \deqn{
    \bar{Y}_i=\frac{\sum_{j=1}^{n_i} Y_{ij}}{n_{i}}, ~~
    \bar{\bm{Y}}=(\bar{Y}_1~ \bar{Y}_2~ \ldots~ \bar{Y}_i)^t, ~~
    T^{\prime}_{k}=\frac{\bm{c}^t_k \bar{\bm{Y}}}{\bm{c}^t_k \bm{c}_k}
  }{
    Ybar_i = (sum from j of Y_ij) / n_i,
    Ybar = (Ybar_1 Ybar_2 ... Ybar_i)^t (i x 1 vector),
    T'_k = c_k^t Ybar / c_k^t c_k
  }
  \deqn{
    T^{\prime}_{\max}=\max(T^{\prime}_1,~ T^{\prime}_2,~ \ldots,~ T^{\prime}_k)
  }{
    T'_max = max(T'_1, T'_2, ..., T'_k)
  }
  
  Consider testing the overall null hypothesis \eqn{H_0: \mu_1=\mu_2=\ldots=\mu_i}, 
  versus alternative hypotheses \eqn{H_1} for response petterns 
  (\eqn{H_1: \mu_1<\mu_2<\ldots<\mu_i,~ \mu_1=\mu_2<\ldots<\mu_i,~ \mu_1<\mu_2<\ldots=\mu_i}).
  The \eqn{P}-value for the probability distribution of \eqn{T^{\prime}_{\max}}{T'_max}
  under the overall null hypothesis is
  \deqn{
    P\mbox{-value}=\Pr(T^{\prime}_{\max}>t^{\prime*}_{\max} \mid H_0)
  }{
    P-value = Pr(T'_max > t*_max | H0)
  }
  \eqn{t^{\prime*}_{\max}}{t*_max} is observed value of statistics.
  This function gives distribution of \eqn{T^{\prime}_{\max}}{T'_max} by using a resampling based procedure.
}
\value{
  \item{cont}{
    suffix of coefficient vector of \eqn{k}th pattern 
    that gives the maximum contrast statistics (colmun number of coefficient matrix).
  }
  \item{pval}{
    \eqn{P}-value of the modified maximum contrast statistics, which is
    simulated by using a resampling based procedure.
  }
  \item{apval}{
    vector of \eqn{P}-values (\eqn{\Pr(T^{\prime}_k>t^{\prime*}_k \mid H_0)}{Pr(T'_max > t*_max | H0})
    for each modified contrast statistics.
  }
  \item{astat}{
    vector of observed values (\eqn{T^{\prime}_k}{T'_k}) for each modified contrast statistics.
  }
  \item{acont}{
    coefficient matrix (copy of argument).
  }
  return value is \code{mmcm.resamp} class object.
}
%\note{
% x list version
% g factor object version
%}
\references{
Sato, Y., Laird, N.M., Nagashima, K., et al. (2009)
A new statistical screening approach for finding pharmacokinetics-related genes in genome-wide studies.
\emph{Pharmacogenomics J.}, in press.
\url{http://www.ncbi.nlm.nih.gov/pubmed/19104505}
}
\author{
Author: Kengo NAGASHIMA and Yasunori SATO

Maintainer: Kengo NAGASHIMA \email{nagasima@josai.ac.jp}
}
\seealso{
  \code{\link{print.mmcm.resamp}}
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
contrast <- cbind(
  c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
)
y <- mmcm.resamp(x, g, contrast, 20000, 5784324)
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
contrast <- cbind(
  c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
)
mmcmtapply <- function(r) {
  mmcm.resamp(
    xx$x[xx$pos==r[1]], xx$g[xx$pos==r[1]],
    contrast, 10000, 5784324+as.numeric(r[1])
  )
}
y <- tapply(xx$pos, xx$pos, mmcmtapply)
yy <- data.frame(
  Pos = as.vector(names(y)),
  Pval = as.vector(sapply(y, "[[", 5)),
  Pattern = as.vector(sapply(y, "[[", 6))
)
yy
}
\keyword{htest}
