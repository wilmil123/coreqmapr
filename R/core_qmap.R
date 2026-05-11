#' Correct predicted core values based on a known top-of-core value
#'
#' `core_qmap` adjusts predicted variables in lake core sediments using a
#' modification of quantile mapping, a climate modelling bias correction
#' technique, using the difference between the predicted value at the top
#' of the core and a known value at the top of the core, usually in water.
#'
#' @details This function returns an object that can make adjustments to the
#' core values using `predict()`, but does not make the adjustments itself.
#' Ensure that the values used in `core_mean`, `core_sd`, `core_top_val`, and
#' `lake_true_val` align between cores.
#'
#' @param core_mean **Numeric vector** of mean values in the cores to be corrected
#' (1 per core).
#' @param core_sd **Numeric vector** of standard deviation values in the cores to
#' be corrected (1 per core).
#' @param core_top_val **Numeric vector** of top-of-core values to be corrected (1
#' per core).
#' @param lake_true_val **Numeric vector** of known values that the top-of-core
#' values should be matched to (1 per core).
#' @param qmap_method **String.**
#' * `"delta"` (default): Quantile mapping method similar to Quantile Delta
#'  Mapping (QDM; Cannon et al., 2015)
#' * `"empirical"`: Quantile mapping method similar to Empirical Quantile
#'  Mapping (EQM)
#' * `"arithmetic`: Not true quantile mapping, but will return transformed
#'  values based on an arithmetic calculation of an optimal value.
#' * `"naive"`: Not true quantile mapping, but will return transformed values
#'  based on a naive scaling of the original values.
#' @param delta_transform **String.** Only used if `qmap_method = "delta"`.
#' * `"multiplicative"` (default): Delta-m should be a multiplicative scalar.
#'  Since the means and standard deviations should be greater than 0, this is
#'  preferable.
#' * `additive`: Delta-m should be an additive scalar.
#' @param interp_method **String.** Any method that can be passed to
#' [imputeTS::na_interpolation()]. Controls how values are imputed for filling
#' out the distribution functions. Default is `"linear"`.
#' @param precision **Integer.** Controls how many decimal places should be used
#' for interpolating distribution functions. Default is `4`.
#' @param conf_alpha **Float.** The alpha value to use for calculating confidence
#' intervals. Default is `0.05`.
#' @param converge_accept **Float.** The acceptance threshold for model convergence.
#' Refers to the minimum error between the modelled prediction and the known
#' value supplied. Default is `0.0001`.
#' @param tempstep **Float.** Controls how quickly the temperature of the searching
#' algorithm decreases, which controls the narrowness of the searching
#' distribution. If not supplied, has different defaults depending on the method.
#' Default for `qmap_method = "delta"` is `0.0001`; default for
#' `qmap_method = "empirical"` is `0.001`.
#' @param iterations **Integer.** Number of maximum iterations to search for optimal
#' weightings, if convergence is not achieved. If not supplied, has different
#' defaults depending on the method. Default for `qmap_method = "delta"` is
#' `10000`; default for `qmap_method = "empirical"` is `100`.
#' @param seed **Integer.** For reproducibility, a seed to pass on to the searching
#' algorithm.
#' @param do_parallel **Boolean.** If `TRUE` (default), samples will be searched
#' for optimal weights in parallel to speed up computation. If `FALSE`, samples
#' will be searched for optimal weights sequentially.
#' @param quiet **Boolean.** If `FALSE` (default), information about weight
#' searching and potential issues will be printed to the console. If `TRUE`,
#' nothing will be printed unless the function encounters a fatal error.
#'
#' @returns An object of class `"coreqm"`, including:
#' \item{`opt_means`}{A vector of optimally transformed means, 1 per core.}
#' \item{`opt_means_ci`}{A data frame including the upper and lower confidence
#'  intervals for the optimally transformed means.}
#' \item{`opt_sds`}{A vector of optimally transformed standard deviations, 1 per
#'  core.}
#' \item{`opt_sds_ci`}{A data frame including the upper and lower confidence
#'  intervals for the optimally transformed standard deviations.}
#' \item{`opt_weights`}{A data frame including the optimal weights for each core.}
#' \item{`internal_data`}{A data frame that includes the original data and
#'  arithmetic and naive transformed values.}
#' \item{`orig_means_ecdf`}{An object of type `"qmrecdf"` for the ECDF of original
#'  core means.}
#' \item{`orig_sds_ecdf`}{An object of type `"qmrecdf"` for the ECDF of original
#'  core standard deviations.}
#' \item{`err_ecdf`}{An object of type `"qmrecdf"` for the ECDF of original
#'  core errors.}
#' \item{`adj_means_ecdf`}{An object of type `"qmrecdf"` for the ECDF of
#'  arithmetically adjusted core means.}
#' \item{`adj_sds_ecdf`}{An object of type `"qmrecdf"` for the ECDF of
#'  arithmetically adjusted core standard deviations.}
#' \item{`naive_means_ecdf`}{An object of type `"qmrecdf"` for the ECDF of naively
#'  adjusted core means.}
#' \item{`naive_sds_ecdf`}{An object of type `"qmrecdf"` for the ECDF of naively
#'  adjusted core standard deviations.}
#' \item{`opt_means_ecdf`}{An object of type `"qmrecdf"` for the ECDF of optimally
#'  adjusted core means.}
#' \item{`opt_sds_ecdf`}{An object of type `"qmrecdf"` for the ECDF of optimally
#'  adjusted core standard deviations.}
#' \item{`weight_error_convergences`}{An object of type `"qdmerrmat"` that contains
#'  information about the convergence of weights as the error function is
#'  minimized.}
#'
#' @export
core_qmap <- function(core_mean,
                      core_sd,
                      core_top_val,
                      lake_true_val,
                      qmap_method = "delta",
                      delta_transform = "multiplicative",
                      interp_method = "linear",
                      precision = 4,
                      conf_alpha = 0.05,
                      converge_accept = 0.0001,
                      tempstep = NULL,
                      iterations = NULL,
                      seed = NULL,
                      do_parallel = TRUE,
                      quiet = FALSE) {
  # set defaults for a few parameters based on which method is being chosen
  if (is.null(tempstep)) {
    if (qmap_method == "delta") {
      tempstep <- 0.0001
    } else if (qmap_method == "empirical") {
      tempstep <- 0.001
    }
  }
  if (is.null(iterations)) {
    if (qmap_method == "delta") {
      iterations <- 10000
    } else if (qmap_method == "empirical") {
      iterations <- 100
    }
  }

  errors_vec <- lake_true_val - core_top_val

  df_to_transform <- data.frame(
    core_mean = core_mean,
    core_sd = core_sd,
    true_val = lake_true_val,
    core_top_val = core_top_val,
    top_err = errors_vec
  )

  df_to_transform$core_top_z <- (df_to_transform$core_top_val - df_to_transform$core_mean) / df_to_transform$core_sd

  adj_mean <- apply(df_to_transform, 1, function(row_number) {
    arithmetic_mean_adjust(row_number["core_mean"], row_number["core_sd"], row_number["core_top_z"], row_number["true_val"])
  })
  df_to_transform$adj_mean <- adj_mean

  adj_sd <- apply(df_to_transform, 1, function(row_number) {
    arithmetic_sd_adjust(row_number["core_mean"], row_number["core_sd"], row_number["core_top_z"], row_number["true_val"])
  })
  df_to_transform$adj_sd <- adj_sd

  naive_mean_sd <- apply(df_to_transform, 1, function(row_number) {
    naive_mean_sd_adjust(row_number["core_mean"], row_number["core_sd"], row_number["top_err"])
  })
  df_to_transform$naive_mean <- naive_mean_sd[, 1]
  df_to_transform$naive_sd <- naive_mean_sd[, 2]

  hist_mean_ecdf <- construct_ecdf_obj(core_mean, precision = precision, method = interp_method)
  hist_sd_ecdf <- construct_ecdf_obj(core_sd, precision = precision, method = interp_method)
  adj_mean_ecdf <- construct_ecdf_obj(adj_mean, precision = precision, method = interp_method)
  adj_sd_ecdf <- construct_ecdf_obj(adj_sd, precision = precision, method = interp_method)
  naive_mean_ecdf <- construct_ecdf_obj(naive_mean_sd[, 1], precision = precision, method = interp_method)
  naive_sd_ecdf <- construct_ecdf_obj(naive_mean_sd[, 2], precision = precision, method = interp_method)
  err_ecdf <- construct_ecdf_obj(errors_vec, precision = precision, method = interp_method)

  if (qmap_method == "arithmetic") {
    opt_transformed_means <- adj_mean
    opt_transformed_sds <- adj_sd
    opt_means_ecdf <- adj_mean_ecdf
    opt_sds_ecdf <- adj_sd_ecdf
    opt_weights <- NA
    weight_error_convergences <- NA
  } else if (qmap_method == "naive") {
    opt_transformed_means <- naive_mean_sd[,1]
    opt_transformed_sds <- naive_mean_sd[,2]
    opt_means_ecdf <- naive_mean_ecdf
    opt_sds_ecdf <- naive_sd_ecdf
    opt_weights <- NA
    weight_error_convergences <- NA
  } else {
    weight_error_convergences <- find_optimal_weight(
      hist_mean_ecdf = hist_mean_ecdf,
      hist_sd_ecdf = hist_sd_ecdf,
      adj_mean_ecdf = adj_mean_ecdf,
      adj_sd_ecdf = adj_sd_ecdf,
      naive_mean_ecdf = naive_mean_ecdf,
      naive_sd_ecdf = naive_sd_ecdf,
      err_ecdf = err_ecdf,
      hist_mean_vals = df_to_transform$core_mean,
      hist_sd_vals = df_to_transform$core_sd,
      adj_mean_vals = df_to_transform$adj_mean,
      adj_sd_vals = df_to_transform$adj_sd,
      naive_mean_vals = df_to_transform$naive_mean,
      naive_sd_vals = df_to_transform$naive_sd,
      err_vals = df_to_transform$top_err,
      z_vals = df_to_transform$core_top_z,
      true_vals = df_to_transform$true_val,
      qmap_method = qmap_method,
      delta_transform = delta_transform,
      tempstep = tempstep,
      iterations = iterations,
      converge_accept = converge_accept,
      quiet = quiet,
      seed = seed,
      do_parallel = do_parallel
    )

    opt_weights <- return_minimum_error(weight_error_convergences)

    if (qmap_method == "empirical") {
      opt_transformed_means <- sapply(1:nrow(df_to_transform), function(row_number) {
        transform_distributions(
          hist_ecdf = hist_mean_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_mean_ecdf,
          naive_ecdf = naive_mean_ecdf,
          hist_val = df_to_transform[row_number, "core_mean"],
          err_val = df_to_transform[row_number, "top_err"],
          adj_val = df_to_transform[row_number, "adj_mean"],
          naive_val = df_to_transform[row_number, "naive_mean"],
          weighting = opt_weights[row_number, "w1"],
          qmap_method = qmap_method,
          delta_transform = delta_transform,
          quiet = quiet
        )
      })

      opt_transformed_sds <- sapply(1:nrow(df_to_transform), function(row_number) {
        transform_distributions(
          hist_ecdf = hist_sd_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_sd_ecdf,
          naive_ecdf = naive_sd_ecdf,
          hist_val = df_to_transform[row_number, "core_sd"],
          err_val = df_to_transform[row_number, "top_err"],
          adj_val = df_to_transform[row_number, "adj_sd"],
          naive_val = df_to_transform[row_number, "naive_sd"],
          weighting = opt_weights[row_number, "w2"],
          qmap_method = qmap_method,
          delta_transform = delta_transform,
          quiet = quiet
        )
      })
    } else if (qmap_method == "delta") {
      opt_transformed_means <- sapply(1:nrow(df_to_transform), function(row_number) {
        transform_distributions(
          hist_ecdf = hist_mean_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_mean_ecdf,
          naive_ecdf = naive_mean_ecdf,
          hist_val = df_to_transform[row_number, "core_mean"],
          err_val = df_to_transform[row_number, "top_err"],
          adj_val = df_to_transform[row_number, "adj_mean"],
          naive_val = df_to_transform[row_number, "naive_mean"],
          weighting = unlist(opt_weights[row_number, c("w1", "w2")], use.names = FALSE),
          qmap_method = qmap_method,
          delta_transform = delta_transform,
          quiet = quiet
        )
      })

      opt_transformed_sds <- sapply(1:nrow(df_to_transform), function(row_number) {
        transform_distributions(
          hist_ecdf = hist_sd_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_sd_ecdf,
          naive_ecdf = naive_sd_ecdf,
          hist_val = df_to_transform[row_number, "core_sd"],
          err_val = df_to_transform[row_number, "top_err"],
          adj_val = df_to_transform[row_number, "adj_sd"],
          naive_val = df_to_transform[row_number, "naive_sd"],
          weighting = unlist(opt_weights[row_number, c("w3", "w4")], use.names = FALSE),
          qmap_method = qmap_method,
          delta_transform = delta_transform,
          quiet = quiet
        )
      })
    }


    opt_means_ecdf <- construct_ecdf_obj(opt_transformed_means,
                                         precision = precision,
                                         method = interp_method)

    opt_sds_ecdf <- construct_ecdf_obj(opt_transformed_sds,
                                       precision = precision,
                                       method = interp_method)
  }

  means_ci <- construct_ci_from_distribution(opt_means_ecdf, opt_transformed_means, conf_alpha)
  sds_ci <- construct_ci_from_distribution(opt_sds_ecdf, opt_transformed_sds, conf_alpha)

  coreqm_returnobj <- list(
    opt_means = opt_transformed_means,
    opt_means_ci = means_ci,
    opt_sds = opt_transformed_sds,
    opt_sds_ci = sds_ci,
    opt_weights = opt_weights,
    internal_data = df_to_transform,
    orig_means_ecdf = hist_mean_ecdf,
    orig_sds_ecdf = hist_sd_ecdf,
    err_ecdf = err_ecdf,
    adj_means_ecdf = adj_mean_ecdf,
    adj_sds_ecdf = adj_sd_ecdf,
    naive_means_ecdf = naive_mean_ecdf,
    naive_sds_ecdf = naive_sd_ecdf,
    opt_means_ecdf = opt_means_ecdf,
    opt_sds_ecdf = opt_sds_ecdf,
    weight_error_convergences = weight_error_convergences
  )
  class(coreqm_returnobj) <- "coreqm"

  return(coreqm_returnobj)
}

