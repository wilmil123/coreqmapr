#' Whole-core predictions using a `coreqm` object
#'
#' This predict method will apply the adjustments from a `coreqm` object to
#' transform downcore values based on the top-of-core adjustments. Variance
#' flattening using 210Pb-inferred dates is also available.
#'
#' @details
#' If either `pb_dates` or `date_col` is (correctly) supplied, the variances
#' of downcore measurements will be flattened according to the number of years
#' contained in each interval in comparison with the first interval. If `pb_dates`
#' is supplied, the flattening will be done according to the full dating profile
#' provided. If `date_col` is supplied, with a properly matching column in
#' `core_vals`, interval widths will be inferred from the modulus of the interval
#' midpoint, and dates will be linearly interpolated up and down to get the date
#' range of the interval.
#' If no variance flattening is desired, leave both `pb_dates` and `date_col`
#' as `NULL`.
#'
#' @param object An object of type `"coreqm"`. See `core_qmap()`, which
#' creates these objects.
#' @param core_vals **A data frame** with at least 3 columns (in order):
#' 1. Column of site/sample IDs.
#' 2. Column of midpoint depths for intervals.
#' 3. Column of original predicted values for each interval.
#' 4. Column of 210Pb dates for each interval (optional).
#' Any columns beyond these first 4 will be ignored.
#' @param pb_dates **A data frame** (optional) that contains a 210Pb dating
#' profile for each core to be adjusted. If supplied, must contain at least 3
#' columns (in order):
#' 1. Column of site/sample IDs to match with `core_vals`.
#' 2. Column of interval top/bottom/midpoint depths.
#' 3. Column of inferred 210Pb ages by interval depth.
#' Any columns beyond the first 3 will be ignored.
#' @param date_col **String.** (optional) If supplying a 4th column to `core_vals`,
#' the name of the column containing 210Pb-inferred dates.
#' @param conf_alpha **Float.** The alpha value to use for calculating confidence
#' intervals. Default is `0.05`.
#' @param quiet **Boolean.** If `FALSE` (default), information about weight
#' searching and potential issues will be printed to the console. If `TRUE`,
#' nothing will be printed unless the function encounters a fatal error.
#' @param clamp_zeroes **Boolean.** Transforming inferred core values in this
#' way can sometimes result in downcore values < 0. If `TRUE` (default), any
#' value below 0 will be clamped to 0. If `FALSE`, no clamping is performed.
#' @param ... Unused.
#'
#' @returns A data frame based on `core_vals` with the following additional
#' columns:
#' \item{transformed_val}{Transformed values of the input measure.}
#' \item{.ci_lower_rng}{Lower range-based C.I.}
#' \item{.ci_upper_rng}{Upper range-based C.I.}
#' \item{.ci_lower_prb}{Lower distribution-based C.I.}
#' \item{.ci_upper_prb}{Upper distribution-based C.I.}
#'
#' @export
predict.coreqm <- function(object,
                           core_vals,
                           pb_dates = NULL,
                           date_col = NULL,
                           conf_alpha = 0.05,
                           quiet = FALSE,
                           clamp_zeroes = TRUE,
                           ...) {
  # within core_vals
  # column 1 is sample ID
  # column 2 is midpoint
  # column 3 is values
  # optional 4th column is for pb date values
  validation_df <- validate_coreqm(object)
  validation_top_ecdf <- construct_ecdf_obj(validation_df$transformed_top_val,
                                            precision = 4,
                                            method = "linear")

  id_col <- core_vals[, 1]
  vals_col <- core_vals[, 3]
  if (!(is.character(id_col) |
        is.factor(core_vals)))
    stop ("The first column must contain a character or factor vector of IDs.")
  if (!(is.numeric(vals_col)))
    stop ("The third column must contain a numeric vector of values")
  site_ids <- unique(id_col)
  transformed_within_site <- lapply(seq_along(site_ids), function(site) {
    current_site <- subset(core_vals, core_vals[, 1] == site_ids[site])
    current_site_dates <- subset(pb_dates, pb_dates[, 1] == site_ids[site])

    within_site_mean <- base::mean(current_site[, 3], na.rm = TRUE)
    within_site_sd <- stats::sd(current_site[, 3], na.rm = TRUE)
    within_core_z <- (current_site[, 3] - within_site_mean) / within_site_sd
    transformed_params_match_id <- lapply(within_site_mean, function(mean) {
      min_diff_to_model_vals <- min(abs(mean - object$internal_data$core_mean))
      if (min_diff_to_model_vals > 0.0001) {
        # warning("The closest site in the model to the current site is more than 0.1 units away. Are you sure this site is properly represented in the data?")
        return(NA)
      }
      id_of_matched_val <- which.min(abs(mean - object$internal_data$core_mean))
      adj_mean <- object$opt_means[id_of_matched_val]
      adj_sd <- object$opt_sds[id_of_matched_val]
      adj_mean_cil <- object$opt_means_ci$ci_lower[id_of_matched_val]
      adj_mean_ciu <- object$opt_means_ci$ci_upper[id_of_matched_val]
      adj_sd_cil <- object$opt_sds_ci$ci_lower[id_of_matched_val]
      adj_sd_ciu <- object$opt_sds_ci$ci_upper[id_of_matched_val]
      params_df <- data.frame(
        adj_mean = adj_mean,
        adj_sd = adj_sd,
        adj_mean_cil = adj_mean_cil,
        adj_mean_ciu = adj_mean_ciu,
        adj_sd_cil = adj_sd_cil,
        adj_sd_ciu = adj_sd_ciu
      )
      return(params_df)
    })
    if (any(is.na(transformed_params_match_id))) {
      early_return_df <- data.frame(
        site_id = current_site[, 1],
        original_val = current_site[, 3],
        transformed_val = NA,
        .ci_lower_rng = NA,
        .ci_upper_rng = NA,
        .ci_lower_prb = NA,
        .ci_upper_prb = NA
      )
      colnames(early_return_df)[1] <- colnames(core_vals)[1]
      colnames(early_return_df)[2] <- colnames(core_vals)[3]
      return(early_return_df)
    } else if (!is.null(pb_dates) &&
               (is.null(current_site_dates) ||
                nrow(current_site_dates) == 0)) {
      early_return_df <- data.frame(
        site_id = current_site[, 1],
        original_val = current_site[, 3],
        transformed_val = NA,
        .ci_lower_rng = NA,
        .ci_upper_rng = NA,
        .ci_lower_prb = NA,
        .ci_upper_prb = NA
      )
      colnames(early_return_df)[1] <- colnames(core_vals)[1]
      colnames(early_return_df)[2] <- colnames(core_vals)[3]
      return(early_return_df)
    }

    if (!is.null(pb_dates) || !is.null(date_col)) {
      if (!is.null(pb_dates)) {
        dateranges <- infer_daterange_with_fulldates(current_site_dates)
      } else if (!is.null(date_col)) {
        dateranges <- infer_daterange_with_datecol(current_site[, 2], current_site[[date_col]])
      }
      corrected_sd_function <- correct_sd(transformed_params_match_id[[1]]$adj_sd, dateranges)
      corrected_sd <- subset(corrected_sd_function,
                             corrected_sd_function$midpoint %in% current_site[, 2])[, 2]
      if (length(corrected_sd) != length(within_core_z))
        stop ("Lengths of adjusted variance function and actual core depths do not match!")

      if (any(corrected_sd < 0)) {
        if (!(quiet)) {
          message_parallel(paste0("[COREQM PREDICTIONS]\n",
                                  "In sample ", site, ":\n",
                                  "Some corrected SD values ended up < 0!"))
        }
      }
    } else {
      corrected_sd <- transformed_params_match_id[[1]]$adj_sd
    }
    untransformed_vals <- (within_core_z * corrected_sd) + transformed_params_match_id[[1]]$adj_mean

    if(any(untransformed_vals < 0)) {
      message_parallel(paste0("[COREQMAP PREDICTIONS]\n",
                              "In sample ", site, ":\n",
                              "Some untransformed values ended up < 0!\n",
                              "This can happen when the values in the original core were already close to 0, or a substantial downwards correction was made."))

      if(clamp_zeroes) {
        if(!(quiet)) {
          message_parallel("Clamping values < 0 to 0.")
        }
        untransformed_vals <- replace(untransformed_vals, untransformed_vals < 0, 0)
      }
    }

    ### ecdf-wise confidence intervals
    vals_ci <- construct_ci_from_distribution(validation_top_ecdf, untransformed_vals, conf_alpha)

    ### range-wise confidence intervals

    if (!is.null(pb_dates) || !is.null(date_col)) {
      if (!is.null(pb_dates)) {
        dateranges <- infer_daterange_with_fulldates(current_site_dates)
      } else if (!is.null(date_col)) {
        dateranges <- infer_daterange_with_datecol(current_site[, 2], current_site[[date_col]])
      }
      corrected_sd_func_cil <- correct_sd(transformed_params_match_id[[1]]$adj_sd_cil, dateranges)
      corrected_sd_cil <- subset(corrected_sd_func_cil,
                                 corrected_sd_func_cil$midpoint %in% current_site[, 2])[, 2]

      corrected_sd_func_ciu <- correct_sd(transformed_params_match_id[[1]]$adj_sd_ciu, dateranges)
      corrected_sd_ciu <- subset(corrected_sd_func_ciu,
                                 corrected_sd_func_ciu$midpoint %in% current_site[, 2])[, 2]
    } else {
      corrected_sd_cil <- transformed_params_match_id[[1]]$adj_sd_cil
      corrected_sd_ciu <- transformed_params_match_id[[1]]$adj_sd_ciu
    }

    vals_ci_rng <- construct_ci_from_sl_range(
      top_z = within_core_z,
      mean_cil = transformed_params_match_id[[1]]$adj_mean_cil,
      mean_ciu = transformed_params_match_id[[1]]$adj_mean_ciu,
      sd_cil = corrected_sd_cil,
      sd_ciu = corrected_sd_ciu
    )

    within_core_df <- data.frame(
      site_id = current_site[, 1],
      original_val = current_site[, 3],
      transformed_val = untransformed_vals,
      .ci_lower_rng = vals_ci_rng$ci_lower,
      .ci_upper_rng = vals_ci_rng$ci_upper,
      .ci_lower_prb = vals_ci$ci_lower,
      .ci_upper_prb = vals_ci$ci_upper
    )
    colnames(within_core_df)[1] <- colnames(core_vals)[1]
    colnames(within_core_df)[2] <- colnames(core_vals)[3]
    return(within_core_df)
  })
  transformed_within_site <- dplyr::bind_rows(transformed_within_site)
  return_df <- dplyr::full_join(core_vals,
                                transformed_within_site,
                                by = c(colnames(core_vals)[1], colnames(core_vals)[3]))
  return(return_df)
}

#' Validate a coreqm object
#'
#' `validate_coreqm` will return the optimal transformed top-of-core values
#' (not downcore) for each sample so that the results can be verified.
#'
#' @param coreqm_obj An object of type `coreqm`.
#' @param interp_method **String.** Any method that can be passed to
#' [imputeTS::na_interpolation()]. Controls how values are imputed for filling
#' out the distribution functions. Default is `"linear"`.
#' @param precision **Integer.** Controls how many decimal places should be used
#' for interpolating distribution functions. Default is `4`.
#' @param conf_alpha **Float.** The alpha value to use for calculating confidence
#' intervals. Default is `0.05`.
#' @param ... Unused.
#'
#' @returns A data frame with the following columns:
#' \item{core_top_val}{The original top-of-core values.}
#' \item{true_val}{The known values.}
#' \item{transformed_top_val}{Modelled top-of-core values.}
#' \item{.ci_lower_rng}{Lower range-based C.I.}
#' \item{.ci_upper_rng}{Upper range-based C.I.}
#' \item{.ci_lower_prb}{Lower distribution-based C.I.}
#' \item{.ci_upper_prb}{Upper distribution-based C.I.}
#'
#' @export
validate_coreqm <- function(coreqm_obj,
                            interp_method = "linear",
                            precision = 4,
                            conf_alpha = 0.05,
                            ...) {
  transformed_top_vals <- (coreqm_obj$internal_data$core_top_z * coreqm_obj$opt_sds) + coreqm_obj$opt_means

  transformed_top_vals_ecdf <- construct_ecdf_obj(transformed_top_vals,
                                                  precision = precision,
                                                  method = interp_method)

  top_vals_ci <- construct_ci_from_distribution(transformed_top_vals_ecdf,
                                                transformed_top_vals,
                                                conf_alpha)

  top_vals_ci_rng <- construct_ci_from_sl_range(
    top_z = coreqm_obj$internal_data$core_top_z,
    mean_cil = coreqm_obj$opt_means_ci$ci_lower,
    mean_ciu = coreqm_obj$opt_means_ci$ci_upper,
    sd_cil = coreqm_obj$opt_sds_ci$ci_lower,
    sd_ciu = coreqm_obj$opt_sds_ci$ci_upper
  )

  pred_df <- data.frame(
    core_top_val = coreqm_obj$internal_data$core_top_val,
    true_val = coreqm_obj$internal_data$true_val,
    transformed_top_val = transformed_top_vals,
    .ci_lower_rng = top_vals_ci_rng$ci_lower,
    .ci_upper_rng = top_vals_ci_rng$ci_upper,
    .ci_lower_prb = top_vals_ci$ci_lower,
    .ci_upper_prb = top_vals_ci$ci_upper
  )
  return(pred_df)
}
