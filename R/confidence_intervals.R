calc_wald_confint <- function(ecdf_obj, prob_i, alpha = 0.05) {
  z_val <- stats::qnorm(1 - (alpha / 2))
  num_obs <- length(which(ecdf_obj$interp_flag == TRUE)) - 2 # subtract 2 for clamping points
  radicand <- (prob_i * (1 - prob_i)) / num_obs
  interv_range <- z_val * sqrt(radicand)
  ci_lower <- prob_i - interv_range
  if (ci_lower < 0) {
    ci_lower <- 0
  }
  ci_upper <- prob_i + interv_range
  if (ci_upper > 1) {
    ci_upper <- 1
  }
  ci_list <- list(ci_lower = ci_lower, ci_upper = ci_upper)
}

construct_ci_from_distribution <- function(ecdf_obj, series, conf_alpha) {
  series_ci <- lapply(series, function(series_item) {
    item_prob <- forward_ecdf_lookup(ecdf_obj, series_item)
    item_prob_ci <- calc_wald_confint(ecdf_obj, item_prob, alpha = conf_alpha)
    item_ci <- data.frame(
      ci_lower = reverse_ecdf_lookup(ecdf_obj, item_prob_ci$ci_lower),
      ci_upper = reverse_ecdf_lookup(ecdf_obj, item_prob_ci$ci_upper)
    )
    return(item_ci)
  })
  series_ci <- dplyr::bind_rows(series_ci)
}

construct_ci_from_sl_range <- function(top_z, mean_cil, mean_ciu, sd_cil, sd_ciu) {
  ci_upper_upper <- (top_z * sd_ciu) + mean_ciu
  ci_upper_lower <- (top_z * sd_ciu) + mean_cil
  ci_lower_upper <- (top_z * sd_cil) + mean_ciu
  ci_lower_lower <- (top_z * sd_cil) + mean_cil

  ci_upper <- sapply(1:length(ci_upper_upper), function(ci_val) {
    max(ci_upper_upper[ci_val],
        ci_upper_lower[ci_val],
        ci_lower_upper[ci_val],
        ci_lower_lower[ci_val])
  })

  ci_lower <- sapply(1:length(ci_upper_upper), function(ci_val) {
    min(ci_upper_upper[ci_val],
        ci_upper_lower[ci_val],
        ci_lower_upper[ci_val],
        ci_lower_lower[ci_val])
  })

  out_df <- data.frame(ci_lower = ci_lower, ci_upper = ci_upper)
}
