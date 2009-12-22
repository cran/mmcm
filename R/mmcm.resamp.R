`mmcm.resamp` <-
function(x, g, contrast, nsample=20000, seed=unclass(Sys.time())) {
  ####################
  # executable check
  ####################
  if (!is.numeric(x)) stop(paste(deparse(substitute(x)), "must be numeric"))
  if (!is.integer(g)) stop(paste(deparse(substitute(g)), "must be integer. ex) g <- as.numeric(factor(g))"))
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
  if (length((1:nrow(contrast))[apply(contrast, 1, sum) != rep(0, nrow(contrast))]) != 0) {
    cat("sum of contrast vector element must be 0\n")
  }
  
  ####################
  # execute mmcm
  ####################
  mmcm_res <- .C(
    "mmcm_rwrap",
    as.double(x),
    as.integer(g),
    as.double(as.vector(t(contrast))),
    as.integer(nsample),
    as.integer(seed),
    as.integer(ncol(contrast)),
    as.integer(nrow(contrast)),
    as.integer(length(x)),
    count=integer(nrow(contrast))
  )
  mmcm_res$acont <- contrast
  mmcm_res$astat <- as.vector(abs(
    (contrast %*% tapply(x, g, mean)) / sqrt(diag(contrast %*% t(contrast))) ))
  mmcm_res$apval <- mmcm_res$count / nsample
  mmcm_res$pval  <- mmcm_res$apval[mmcm_res$astat==max(mmcm_res$astat)]
  mmcm_res$cont  <- (1:length(mmcm_res$apval))[mmcm_res$astat==max(mmcm_res$astat)]
  mmcm_res[1:8]  <- NULL

  ####################
  # return result
  ####################
  class(mmcm_res) <- "mmcm.res"
  return(mmcm_res)

}

