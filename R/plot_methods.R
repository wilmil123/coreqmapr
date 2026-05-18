#' Plot Interpolated ECDFs
#'
#' Plot method for class `"qmrecdf"`. Capable of plotting the interpolated
#' ECDF along with a confidence interval.
#'
#' @param x An object of class `"qmrecdf"`.
#' @param ci **Boolean.** If `TRUE` (default), plot C.I. around ECDF values.
#' @param conf_alpha **Float.** The alpha value to use for calculating confidence
#' intervals. Default is `0.05`.
#' @param col **String.** Colour to plot confidence interval lines.
#' @param ... Further arguments to pass to `base::plot()`.
#' @export
plot.qmrecdf <- function(x,
                         ci = TRUE,
                         conf_alpha = 0.05,
                         col = "blue",
                         ...) {
  if (ci) {
    x$p_lower_ci <- apply(x, 1, function(row) {
      return(calc_wald_confint(x, row["prob"], alpha = conf_alpha)[[1]])
    })
    x$p_upper_ci <- apply(x, 1, function(row) {
      return(calc_wald_confint(x, row["prob"], alpha = conf_alpha)[[2]])
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

#' Plot Convergence-Error values for Simulated Annealing Results
#'
#' Plot method for objects of class `"qmrerrmat"`.
#'
#' @param x An object of class `"qmrerrmat"`.
#' @param best_col **String.** Colour for the best values.
#' @param try_col **String.** Colour for the tried (neighbour) values.
#' @param ... Further arguments to be passed to `base::plot()`.
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
