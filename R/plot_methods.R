#' @export
plot.qmrecdf <- function(x,
                         ci = TRUE,
                         col = "blue",
                         ...) {
  if (ci) {
    x$p_lower_ci <- apply(x, 1, function(row) {
      return(calc_wald_confint(x, row["prob"])[[1]])
    })
    x$p_upper_ci <- apply(x, 1, function(row) {
      return(calc_wald_confint(x, row["prob"])[[2]])
    })
  }

  plot(x$value,
       x$prob,
       type = "l",
       xlab = "value",
       ylab = "prob",
       ...)
  if (ci) {
    graphics::lines(x$value, x$p_lower_ci, type = "l", col = col)
    graphics::lines(x$value, x$p_upper_ci, type = "l", col = col)
  }
}

#' @export
plot.qmrerrmat <- function(x,
                           best_col = "red",
                           try_col = "blue",
                           ...) {
  plot(
    x$iteration,
    log10(sqrt(x$best_fit)),
    type = "l",
    col = best_col,
    xlab = "iteration",
    ylab = "log10 RSerr",
    ...
  )
  graphics::lines(x$iteration, log10(sqrt(x$fit)), type = "l", col = try_col)
}
