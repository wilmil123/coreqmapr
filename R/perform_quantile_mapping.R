transform_distributions <- function(hist_ecdf,
                                    err_ecdf,
                                    adj_ecdf,
                                    naive_ecdf,
                                    hist_val,
                                    err_val,
                                    adj_val,
                                    naive_val,
                                    weighting,
                                    qmap_method,
                                    delta_transform,
                                    quiet) {
  switch(qmap_method,
         "delta" = {
           if (length(weighting) != 2)
             stop ("Wrong weighting length supplied.")
           err_quantile <- forward_ecdf_lookup(err_ecdf, err_val)
           hist_quantile <- forward_ecdf_lookup(hist_ecdf, hist_val)
           naive_quantile <- forward_ecdf_lookup(naive_ecdf, naive_val)
           arith_quantile <- forward_ecdf_lookup(adj_ecdf, adj_val)
           adj_quantile <- (hist_quantile * (1 - weighting[1])) + (err_quantile * weighting[1])

           initial_proj <- reverse_ecdf_lookup(adj_ecdf, adj_quantile)

           if (delta_transform == "multiplicative") {
             delta_m <- naive_val / reverse_ecdf_lookup(adj_ecdf, naive_quantile)

             if (!(quiet) && delta_m < 0) {
               warning("DELTAm ended up < 0!")
               message_parallel(
                 paste0(
                   "[TRANSFORM DISTRIBUTIONS]\n",
                   "DELTAm ",
                   delta_m,
                   "\n",
                   "ARITH VAL ",
                   adj_val,
                   "\n",
                   "NAIVE VAL ",
                   naive_val,
                   "\n\n"
                 )
               )
             }

             corrected_proj <- initial_proj * (delta_m ** weighting[2])
           } else if (delta_transform == "additive") {
             ### DON'T WORRY ABOUT THIS FOR NOW
             delta_m <- naive_val - reverse_ecdf_lookup(adj_ecdf, naive_quantile)
             corrected_proj <- initial_proj + (delta_m * weighting[2])
           } else {
             stop("Unknown delta transform type.")
           }

           return(corrected_proj)
         },
         "empirical" = {
           if (length(weighting) != 1)
             stop ("Wrong weighting length supplied.")
           err_quantile <- forward_ecdf_lookup(err_ecdf, err_val)
           hist_quantile <- forward_ecdf_lookup(hist_ecdf, hist_val)
           adj_quantile <- (hist_quantile * (1 - weighting)) + (err_quantile * weighting)
           corrected_proj <- reverse_ecdf_lookup(adj_ecdf, adj_quantile)
           return(corrected_proj)
         },
         stop("Something went wrong :("))
}

find_optimal_weight <- function(hist_mean_ecdf,
                                hist_sd_ecdf,
                                adj_mean_ecdf,
                                adj_sd_ecdf,
                                naive_mean_ecdf,
                                naive_sd_ecdf,
                                err_ecdf,
                                hist_mean_vals,
                                hist_sd_vals,
                                adj_mean_vals,
                                adj_sd_vals,
                                naive_mean_vals,
                                naive_sd_vals,
                                err_vals,
                                z_vals,
                                true_vals,
                                qmap_method,
                                delta_transform,
                                tempstep,
                                iterations,
                                converge_accept,
                                quiet,
                                seed,
                                do_parallel) {
  stopifnot(length(hist_mean_vals) == length(hist_sd_vals))

  ### there has to be a better way to do this
  ### consider refactoring

  simulate_annealing <- function(sample) {
    # simulated annealing to find best weights
    # https://codereview.stackexchange.com/questions/84688/simulated-annealing-in-r

    # initialization of tracking vectors
    out_iterations <- vector(length = iterations)
    out_w1 <- vector(length = iterations)
    out_w2 <- vector(length = iterations)
    out_fits <- vector(length = iterations)
    out_bestfits <- vector(length = iterations)
    out_temps <- vector(length = iterations)

    if (qmap_method == "empirical") {
      initial_state <- c(0, 0) # w1, w2
    } else if (qmap_method == "delta") {
      initial_state <- c(0, 0, 0, 0) # w1, w2, w3, w4
      out_w3 <- out_w4 <- vector(length = iterations)
    }

    best_state <- current_state <- neighbour_state <- initial_state

    if (qmap_method == "empirical") {
      transformed_mean <- transform_distributions(
        hist_ecdf = hist_mean_ecdf,
        err_ecdf = err_ecdf,
        adj_ecdf = adj_mean_ecdf,
        naive_ecdf =  naive_mean_ecdf,
        hist_val = hist_mean_vals[sample],
        err_val = err_vals[sample],
        adj_val = adj_mean_vals[sample],
        naive_val = naive_mean_vals[sample],
        weighting = initial_state[1],
        qmap_method = qmap_method,
        delta_transform = delta_transform,
        quiet = quiet
      )

      transformed_sd <- transform_distributions(
        hist_ecdf = hist_sd_ecdf,
        err_ecdf = err_ecdf,
        adj_ecdf = adj_sd_ecdf,
        naive_ecdf = naive_mean_ecdf,
        hist_val = hist_sd_vals[sample],
        err_val = err_vals[sample],
        adj_val = adj_sd_vals[sample],
        naive_val = naive_sd_vals[sample],
        weighting = initial_state[2],
        qmap_method = qmap_method,
        delta_transform = delta_transform,
        quiet = quiet
      )
    } else if (qmap_method == "delta") {
      transformed_mean <- transform_distributions(
        hist_ecdf = hist_mean_ecdf,
        err_ecdf = err_ecdf,
        adj_ecdf = adj_mean_ecdf,
        naive_ecdf = naive_mean_ecdf,
        hist_val = hist_mean_vals[sample],
        err_val = err_vals[sample],
        adj_val = adj_mean_vals[sample],
        naive_val = naive_mean_vals[sample],
        weighting = initial_state[1:2],
        qmap_method = qmap_method,
        delta_transform = delta_transform,
        quiet = quiet
      )

      transformed_sd <- transform_distributions(
        hist_ecdf = hist_sd_ecdf,
        err_ecdf = err_ecdf,
        adj_ecdf = adj_sd_ecdf,
        naive_ecdf = naive_sd_ecdf,
        hist_val = hist_sd_vals[sample],
        err_val = err_vals[sample],
        adj_val = adj_sd_vals[sample],
        naive_val = naive_sd_vals[sample],
        weighting = initial_state[3:4],
        qmap_method = qmap_method,
        delta_transform = delta_transform,
        quiet = quiet
      )
    }


    untransformed_val <- (z_vals[sample] * transformed_sd) + transformed_mean
    sq_err <- (untransformed_val - true_vals[sample]) ** 2

    best_value <- current_value <- neighbour_value <- sq_err

    if (!is.null(seed))
      set.seed(seed)

    for (cur_iter in 1:iterations) {
      #############
      ### DEBUG ###
      #############
      if (!(quiet)) {
        message_parallel(
          paste0(
            "[SIMULATED ANNEALING]\n",
            "Sample ",
            sample,
            " of ",
            length(hist_mean_vals),
            "\n",
            "Iteration ",
            cur_iter
          )
        )
      }
      #############

      temperature <- (1 - tempstep) ** cur_iter

      # https://stats.stackexchange.com/questions/316086/distribution-that-has-a-range-from-0-to-1-and-with-peak-between-them
      # this means the distribution will become narrower as the temperature increases
      narrowness <- 1 + (1 / temperature)

      neighbour_state <- purrr::map_dbl(current_state, function(mean) {
        if (mean == 0) {
          mean <- mean + 0.001
        }
        if (mean < 0.5) {
          alpha <- narrowness
          beta <- ((-alpha * mean) + alpha) / mean
        } else {
          beta <- narrowness
          alpha <- (-beta * mean) / (mean - 1)
        }

        nbr <- stats::rbeta(
          n = 1,
          shape1 = alpha + 1,
          shape2 = beta + 1,
          ncp = 0
        )
      })

      # compute outcome

      if (qmap_method == "empirical") {
        transformed_mean <- transform_distributions(
          hist_ecdf = hist_mean_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_mean_ecdf,
          naive_ecdf = naive_mean_ecdf,
          hist_val = hist_mean_vals[sample],
          err_val = err_vals[sample],
          adj_val = adj_mean_vals[sample],
          naive_val = naive_mean_vals[sample],
          weighting = neighbour_state[1],
          qmap_method = qmap_method,
          delta_transform = delta_transform
        )

        transformed_sd <- transform_distributions(
          hist_ecdf = hist_sd_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_sd_ecdf,
          naive_ecdf = naive_sd_ecdf,
          hist_val = hist_sd_vals[sample],
          err_val = err_vals[sample],
          adj_val = adj_sd_vals[sample],
          naive_val = naive_sd_vals[sample],
          weighting = neighbour_state[2],
          qmap_method = qmap_method,
          delta_transform = delta_transform
        )
      } else if (qmap_method == "delta") {
        transformed_mean <- transform_distributions(
          hist_ecdf = hist_mean_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_mean_ecdf,
          naive_ecdf = naive_mean_ecdf,
          hist_val = hist_mean_vals[sample],
          err_val = err_vals[sample],
          adj_val = adj_mean_vals[sample],
          naive_val = naive_mean_vals[sample],
          weighting = neighbour_state[1:2],
          qmap_method = qmap_method,
          delta_transform = delta_transform,
          quiet = quiet
        )

        transformed_sd <- transform_distributions(
          hist_ecdf = hist_sd_ecdf,
          err_ecdf = err_ecdf,
          adj_ecdf = adj_sd_ecdf,
          naive_ecdf = naive_sd_ecdf,
          hist_val = hist_sd_vals[sample],
          err_val = err_vals[sample],
          adj_val = adj_sd_vals[sample],
          naive_val = naive_sd_vals[sample],
          weighting = neighbour_state[3:4],
          qmap_method = qmap_method,
          delta_transform = delta_transform,
          quiet = quiet
        )
      }

      untransformed_val <- (z_vals[sample] * transformed_sd) + transformed_mean
      neighbour_value <- (untransformed_val - true_vals[sample]) ** 2

      # update current state
      if (neighbour_value < current_value ||
          stats::runif(n = 1, min = 0, max = 1) < exp(-(neighbour_value - current_value) / temperature)) {
        current_state <- neighbour_state
        current_value <- neighbour_value
      }

      # update best state
      if (neighbour_value < best_value) {
        best_state <- neighbour_state
        best_value <- neighbour_value
      }

      # track all the values
      out_iterations[cur_iter] <- cur_iter
      out_w1[cur_iter] <- current_state[1]
      out_w2[cur_iter] <- current_state[2]

      if (qmap_method == "delta") {
        out_w3[cur_iter] <- current_state[3]
        out_w4[cur_iter] <- current_state[4]
      }

      out_fits[cur_iter] <- neighbour_value
      out_bestfits[cur_iter] <- best_value
      out_temps[cur_iter] <- temperature

      # if our acceptance threshold is reached, end here
      if (best_value <= converge_accept) {
        out_iterations <- out_iterations[1:cur_iter]
        out_w1 <- out_w1[1:cur_iter]
        out_w2 <- out_w2[1:cur_iter]

        if (qmap_method == "delta") {
          out_w3 <- out_w3[1:cur_iter]
          out_w4 <- out_w4[1:cur_iter]
        }

        out_fits <- out_fits[1:cur_iter]
        out_bestfits <- out_bestfits[1:cur_iter]
        out_temps <- out_temps[1:cur_iter]
        # DEBUG
        if (!(quiet)) {
          message_parallel(
            paste0(
              "[SIMULATED ANNEALING]\n",
              "Sample ",
              sample,
              " converged after ",
              cur_iter,
              " iterations."
            )
          )
        }
        # END DEBUG
        break
      }
    }
    if (qmap_method == "empirical") {
      err_mat <- data.frame(
        iteration = out_iterations,
        w1 = out_w1,
        w2 = out_w2,
        fit = out_fits,
        best_fit = out_bestfits,
        temperature = out_temps
      )
    } else if (qmap_method == "delta") {
      err_mat <- data.frame(
        iteration = out_iterations,
        w1 = out_w1,
        w2 = out_w2,
        w3 = out_w3,
        w4 = out_w4,
        fit = out_fits,
        best_fit = out_bestfits,
        temperature = out_temps
      )
    }

    class(err_mat) <- c("qmrerrmat", class(err_mat))
    return(err_mat)
  }

  if (do_parallel) {
    err_list <- parallel::mclapply(X = 1:length(hist_mean_vals), FUN = simulate_annealing)
  } else {
    err_list <- lapply(1:length(hist_mean_vals), simulate_annealing)
  }

  return(err_list)
}

return_minimum_error <- function(err_list) {
  opt_weights <- lapply(seq_along(err_list), function(id) {
    sample <- err_list[[id]]
    id_which_min_err <- which.min(sample$best_fit)
    opt_weights_within_sample <- dplyr::select(sample[id_which_min_err, ], tidyselect::starts_with("w"))
    opt_weights_within_sample$min_err <- sample[id_which_min_err, "best_fit"]
    return(opt_weights_within_sample)
  })
  opt_weights <- dplyr::bind_rows(opt_weights)
  return(opt_weights)
}
