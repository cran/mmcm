`mmcm.resamp` <-
function(x, g, contrast, nsample=20000, seed=unclass(Sys.time())) {
  ####################
  # executable check
  ####################
  if (!is.numeric(x)) stop(paste(deparse(substitute(x)), "must be numeric"))
  if (!is.integer(g)) stop(paste(deparse(substitute(g)), "must be integer"))
  if (!is.numeric(contrast)) stop(paste(deparse(substitute(contrast)), "must be numeric"))
  if (!is.numeric(nsample)) stop("nsample must be integer")
  if (!is.numeric(seed)) stop("seed must be integer")
  if (length(x) <= 0) {
    stop(paste(deparse(substitute(x)), " length (", length(x), ") must be > 0", sep=""))
  }
  if (nsample <= 0) stop("nsample must be > 0")
  if (seed <= -1 * 2^31 | seed >= 2^31) stop("seed is a 32-bit integer (-1 * 2^31 + 1 <= seed <= 2^31 - 1)")
  if (length(x) != length(g)) {
    stop(paste(
      deparse(substitute(x)), " length (", length(x), ") is not equal ",
      deparse(substitute(g)), " length (", length(g), ")", sep=""))
  }
  if (length(unique(g)) != nrow(contrast)) {
    stop(paste("dimension of contrast vector is small\n  No. of group =", length(unique(g))))
  }
  if (length((1:ncol(contrast))[apply(contrast, 2, sum) != rep(0, ncol(contrast))]) != 0) {
    cat("sum of contrast vector element must be 0\n")
  }

  ####################
  # execute mmcm
  ####################
  mmcm_res <- .C(
    "mmcm_rwrap",
    as.double(x),
    as.integer(g),
    as.double(as.vector(contrast)),
    as.integer(nsample),
    as.integer(seed),
    as.integer(nrow(contrast)),
    as.integer(ncol(contrast)),
    as.integer(length(x)),
    count=integer(ncol(contrast))
  )
  mmcm_res$acont <- contrast
  mmcm_res$astat <- as.vector(abs(
    (t(contrast) %*% tapply(x, g, mean)) / sqrt(diag(t(contrast) %*% contrast)) ))
  mmcm_res$apval <- mmcm_res$count / nsample
  mmcm_res$pval  <- mmcm_res$apval[mmcm_res$astat==max(mmcm_res$astat)]
  mmcm_res$cont  <- (1:length(mmcm_res$apval))[mmcm_res$astat==max(mmcm_res$astat)]
  mmcm_res[1:8]  <- NULL

  ####################
  # return result
  ####################
  class(mmcm_res) <- "mmcm.resamp"
  return(mmcm_res)

}

