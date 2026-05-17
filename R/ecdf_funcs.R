construct_ecdf_obj <- function(in_vec, precision, method) {
  ecdf_obj <- structure(interpolate_ecdf(expand_ecdf(calc_custom_ecdf(in_vec), precision = precision), method = method),
                        class = c("qmrecdf", "data.frame"))
}

calc_custom_ecdf <- function(in_vec) {
  vec_inorder <- sort(in_vec)
  distributions <- vapply(seq_along(vec_inorder), function(index) {
    return(sum(vec_inorder <= vec_inorder[index]) / length(in_vec))
  }, double(1))
  out_df <- data.frame(value = vec_inorder, prob = distributions)
}

expand_ecdf <- function(ecdf_obj, precision = 4) {
  # set last point = ~0.99
  # approximated based on n values by 1 - (1/(n * precision))
  # for e.g. 30 samples, precision 4, yields 1 - (1/120)
  # ~ 0.99167
  ecdf_obj$prob[length(ecdf_obj$prob)] <- 1 - (1 / (length(ecdf_obj$value) * precision))
  # clamp 0 to min - 0.1*range
  ecdf_obj <- rbind(c(
    value = min(ecdf_obj$value) - 0.1 * (max(ecdf_obj$value) - min(ecdf_obj$value)),
    prob = 0
  ), ecdf_obj)
  # clamp 1 to max + 0.1*range
  ecdf_obj <- rbind(ecdf_obj, c(
    value = max(ecdf_obj$value) + 0.1 * (max(ecdf_obj$value) - min(ecdf_obj$value)),
    prob = 1
  ))
  ecdf_obj$prob <- round(ecdf_obj$prob, precision)

  expanded_ecdf <- tidyr::complete(ecdf_obj, prob = round(seq(0, 1, by = 1 /
                                                                (10 ** precision)), precision))
}

interpolate_ecdf <- function(exp_ecdf_obj, method = "linear") {
  interp_ecdf <- exp_ecdf_obj
  interp_ecdf$value <- imputeTS::na_interpolation(interp_ecdf$value, option = method)
  interp_ecdf$interp_flag <- ifelse(interp_ecdf$value %in% exp_ecdf_obj$value, TRUE, FALSE)
}

forward_ecdf_lookup <- function(ecdf_obj, value) {
  lookup_val <- ecdf_obj$prob[which.min(abs(ecdf_obj$value - value))]
}

reverse_ecdf_lookup <- function(ecdf_obj, prob) {
  lookup_val <- ecdf_obj$value[which.min(abs(ecdf_obj$prob - prob))]
}
