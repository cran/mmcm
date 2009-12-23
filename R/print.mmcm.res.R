`print.mmcm.res` <-
function(x, ...) {
  ####################
  # print result
  ####################
  if (length(unique(x$astat))!=length(x$astat)) {
    cat(sprintf("\n\tModified maximum contrast method\n\n"))
    cat("more than 2 contrast vectors was selected ...\n\n")
    cat(
      sprintf("selected contrast = "),
      sprintf("%d ", x$cont),
      sprintf("\nP-values = "),
      sprintf("%f ", x$pval),
      sprintf("\b\n\n"), sep=""
    )
  } else {
    cat(
      sprintf("\n\t%s\n\n", x$method),
      sprintf("observed statistics = %.3f\n", x$astat[x$cont]),
      sprintf("contrast = %d\t(", x$cont),
      sprintf("%s ", x$acont[,x$cont]),
      sprintf("\b)\n"),
      sprintf("P-value = %f\n\n", x$pval), sep=""
    )
  }
  invisible(0)
}

