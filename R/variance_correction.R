infer_daterange_with_datecol <- function(midpoint_col, age_col) {
  interval_widths <- ifelse((midpoint_col * 2) %% 1 == 0, 1, (midpoint_col * 2) %% 1)
  interval_tops <- midpoint_col - (interval_widths / 2)
  interval_bottoms <- midpoint_col + (interval_widths / 2)
  expanded_intervals <- unique(sort(c(interval_tops, midpoint_col, interval_bottoms)))
  min_interval_dist <- min(diff(expanded_intervals))
  expanded_intervals_tointerp <- seq(min(expanded_intervals),
                                     max(expanded_intervals),
                                     by = min_interval_dist)

  core_series <- data.frame(midpoint = midpoint_col, age = age_col)

  expanded_midpoints <- tidyr::expand(core_series, interp_midpoint = expanded_intervals_tointerp)
  join_cols <- c("interp_midpoint" = "midpoint")
  expanded_data <- dplyr::left_join(expanded_midpoints, core_series, by = join_cols)
  expanded_data$interp_year <- imputeTS::na_interpolation(expanded_data$age, option = "linear")
  expanded_data <- subset(expanded_data, expanded_data$interp_midpoint %in% expanded_intervals)

  expanded_data$year_prev <- c(expanded_data$interp_year[-1], expanded_data$interp_year[length(expanded_data$interp_year)])
  expanded_data$year_next <- c(expanded_data$interp_year[1], expanded_data$interp_year[-length(expanded_data$interp_year)])
  expanded_data$age_diffs <- expanded_data$year_next - expanded_data$year_prev
  expanded_data_subset <- subset(expanded_data,
                                 expanded_data$interp_midpoint %in% midpoint_col)
  age_diffs_df <- data.frame(midpoint = expanded_data_subset$interp_midpoint,
                             age_diff = expanded_data_subset$age_diffs)
}

infer_daterange_with_fulldates <- function(pbdates) {
  # column 1 should be ids
  # column 2 should be interval depths
  # column 3 should be ages
  interval_width <- pbdates[2, 2] - pbdates[1, 2]
  full_intervals <- subset(pbdates, pbdates[, 2] %% (2 * interval_width) == 0)
  half_intervals <- subset(pbdates, pbdates[, 2] %% (2 * interval_width) != 0)
  age_diffs <- diff(full_intervals[, 3])

  if (pbdates[nrow(pbdates), 2] %% (2 * interval_width) != 0) {
    # duplicate the last value if the lengths are wrong
    age_diffs <- c(age_diffs, age_diffs[length(age_diffs)])
  }

  age_diffs_df <- data.frame(midpoint = half_intervals[, 2], age_diff = age_diffs)
}

correct_sd <- function(sdev, age_diffs) {
  variance <- sdev ** 2
  scalefactor <- age_diffs$age_diff[1] / age_diffs$age_diff
  scaled_variance <- scalefactor * variance
  corrected_sd <- sqrt(scaled_variance)
  matched_df <- data.frame(midpoint = age_diffs$midpoint, corrected_sd = corrected_sd)
  print(matched_df)
}
