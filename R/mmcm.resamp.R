`mmcm.resamp` <-
function(x, g, contrast, alternative = c("two.sided", "less", "greater"),
  nsample=20000, abseps=0.001, seed=NULL) {
  
  ####################
  # executable check
  ####################
  
  alternative <- match.arg(alternative)
  
  DNAMEX <- deparse(substitute(x))
  DNAMEG <- deparse(substitute(g))
  DNAMEC <- deparse(substitute(contrast))
  DNAME  <- paste("'", DNAMEX, "' by group '", DNAMEG,
                  "' with contrast coefficient matrix '",
                  DNAMEC, "'", sep="")
  
  if (!is.numeric(x)) {
    stop(paste(DNAMEX, "must be numeric"))
  }
  if (!is.numeric(g)) {
    stop(paste(DNAMEG, "must be numeric"))
  }
  if (!is.matrix(contrast)) {
    stop(paste(DNAMEC, "must be a matrix"))
  }
  
  x <- x[is.finite(x)]
  g <- g[is.finite(g)]
  
  if (length(x) < 1L) {
    stop(paste("not enough (finite) ", DNAMEX, "observations"))
  }
  
  if (length(x) != length(g)) {
    stop(paste(DNAME, "and", DNAMEG, "must have the same length"))
  }
  
  if (length(unique(g)) != ncol(contrast)) {
    stop(paste("nrow(", DNAMEC, ") and length(unique(", DNAMEG,
               ")) must have the same length", sep=""))
  }
  
  if (length((1:nrow(contrast))[apply(contrast, 1, sum) != rep(0, nrow(contrast))]) != 0) {
    stop("sum of contrast vector element must be 0\n")
  }
  
  ####################
  # execute mmcm
  ####################
  
  METHOD <- "Permuted modified maximum contrast method"
  
  p          <- length(unique(g))
  m          <- nrow(contrast)
  df         <- length(g) - p
  pooled     <- (tapply(x, g, length) - 1) * tapply(x, g, var)
  pooled     <- t(rep(1, p)) %*% pooled / df
  D          <- diag(1 / as.vector(tapply(x, g, length)))
  CDC        <- contrast %*% D %*% t(contrast)
  
  CtC        <- contrast %*% t(contrast)
  CtCMATRIX  <- matrix(rep(1, m * m), ncol=m) * diag(CtC)
  Rs         <- CDC / (sqrt(CtCMATRIX) * sqrt(t(CtCMATRIX)))
  
  STATISTICS <- switch(
    alternative,
    less      = (contrast %*% tapply(x, g, mean) / sqrt(diag(CtC))),
    greater   = (contrast %*% tapply(x, g, mean) / sqrt(diag(CtC))),
    two.sided = abs(contrast %*% tapply(x, g, mean) / sqrt(diag(CtC)))
  )
  STATISTIC  <- switch(
    alternative,
    less      = min(STATISTICS),
    greater   = max(STATISTICS),
    two.sided = max(STATISTICS)
  )
  IMAXCONT <- (1:m)[STATISTICS==STATISTIC]
  NMAXCONT <- contrast[IMAXCONT,]
  
  nalternative <- switch(
    alternative,
    less      = 1,
    greater   = 2,
    two.sided = 3
  )
  if (!is.null(seed)) {
    set.seed(seed)
  }
  RESAMP <- .C(
    "mmcm_rwrap",
    as.double(x),
    as.double(g),
    as.double(as.vector(t(contrast))),
    as.integer(nsample),
    as.integer(ncol(contrast)),
    as.integer(nrow(contrast)),
    as.integer(length(x)),
    as.double(abseps),
    as.integer(nalternative),
    pval=double(1),
    error=double(1)
  )
  
  if (RESAMP$error < abseps) {
    msg <- "Normal Completion"
  } else {
    msg <- warning("Completion with error > abseps")
  }
  
  PVAL  <- RESAMP$pval
  ERROR <- RESAMP$error
  MSG   <- msg
  
  if (length((1:m)[STATISTICS==STATISTIC]) != 1) {
    MAXCONT <- warning("More than 2 contrast coefficient vectors were selected")
  } else {
    MAXCONT <- "("
    for(i in 1:p) {
      if (i==p) {
        MAXCONT <- paste(MAXCONT, NMAXCONT[i], ")", sep="")
      } else {
        MAXCONT <- paste(MAXCONT, NMAXCONT[i], ", ", sep="")
      }
    }
  }
  
  names(STATISTIC) <- "Permuted modified maximum contrast statistic"
  names(IMAXCONT)  <- "index"
  names(MAXCONT)   <- "Maximum contrast coefficient vector"
  names(ERROR)     <- "Estimated absolute error of P-value"
  names(MSG)       <- "Status messages of P-value calculation"
  
  RVAL <- structure(list(
    statistic      = STATISTIC,
    parameter      = NULL,
    p.value        = as.numeric(PVAL), 
    alternative    = alternative,
    method         = METHOD, 
    data.name      = DNAME,
    contrast       = MAXCONT,
    contrast.index = IMAXCONT,
    error          = ERROR,
    msg            = MSG),
    class          = "mmcm"
  )
  return(RVAL)
  
}
