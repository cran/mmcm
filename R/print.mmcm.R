`print.mmcm` <-
function(x, digits = getOption("digits"), ...) {

  ####################
  # print result
  ####################
  if (x$contrast=="More than 2 contrast coefficient vectors were selected") {
    class(x) <- "htest"
    print(x)
    cat("More than 2 contrast coefficient vectors were selected\n\n")
  } else {
    class(x) <- "htest"
    print(x)
    msg <- paste(
      names(x$contrast      ), " = ", x$contrast      , ", ",
      names(x$contrast.index), " = ", x$contrast.index, "\n\n",
      names(x$error         ), " = ", format.pval(x$error, digits = digits), "\n",
      names(x$msg           ), " = ", x$msg           , "\n\n",
      sep = ""
    )
    cat(msg)
  }
  invisible(0)
}

