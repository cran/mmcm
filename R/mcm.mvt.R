`mcm.mvt` <-
function(x, g, contrast, algorithm = GenzBretz()) {
  ####################
  # executable check
  ####################
  if (!is.numeric(x)) stop(paste(deparse(substitute(x)), "must be numeric"))
  if (!is.integer(g)) stop(paste(deparse(substitute(g)), "must be integer. ex) g <- as.numeric(factor(g))"))
  if (!is.numeric(contrast)) stop(paste(deparse(substitute(contrast)), "must be numeric"))
  if (length(x) != length(g)) {
    stop(paste(
      deparse(substitute(x)), " length (", length(x), ") is not equal ",
      deparse(substitute(g)), " length (", length(g), ")", sep=""))
  }
  if (length(unique(g)) != nrow(contrast)) {
    stop(paste("dimension of contrast vector is small\n  No. of group =", length(unique(g))))
  }
  if (length((1:nrow(contrast))[apply(contrast, 1, sum) != rep(0, nrow(contrast))]) != 0) {
    cat("sum of contrast vector element must be 0\n")
  }

  ####################
  # execute mmcm
  ####################
  mmcm_res <- NULL

  p  <- length(unique(g))
  m  <- nrow(contrast)
  df <- length(g) - p
  pooled <- (tapply(x, g, length)-1) * tapply(x, g, var)
  pooled <- t(rep(1, p)) %*% pooled / df
  D <- diag(1 / as.vector(tapply(x, g, length)))

  CtC <- contrast %*% t(contrast)
  CDC <- contrast %*% D %*% t(contrast)

  B3 <- diag(1/sqrt(diag(CDC)))
  R  <- B3 %*% CDC %*% B3

  mmcm_res$acont <- contrast

  # sample statistics
  mmcm_res$astat <- abs( contrast %*% tapply(x, g, mean) / sqrt(diag(CDC)) ) / (rep(1, m) %*% sqrt(pooled))
  tmax <- max(mmcm_res$astat)

  # p-value
  mmcm_res$apval <- NULL
  mmcm_res$pval  <- 1 - pmvt(lower = rep(-tmax, p), upper = rep(tmax, p), df = df, sigma = R, algorithm = algorithm)

  # detect contrast
  mmcm_res$cont <- (1:m)[mmcm_res$astat==tmax]

  mmcm_res$method <- "the maximum contrast method"

  ####################
  # return result
  ####################
  class(mmcm_res) <- "mmcm.res"
  return(mmcm_res)

}

