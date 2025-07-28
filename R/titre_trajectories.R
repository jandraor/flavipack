simulate_titre_trajectory <- function(sampling_times, serostatus,
                                      treatment_group,
                                      subject_id,
                                      vac_times     = NULL,
                                      symp_inf_time = NULL)
{
  n_samples     <- length(sampling_times)
  enrolment_day <- min(sampling_times)
  dropout_day   <- max(sampling_times)

  max_inf <- ifelse(serostatus == 0, 4, 3)

  set.seed(1652)

  infection_times <- simulate_infection_times(
    lambda    = 0.3 / 365,
    stop_time = dropout_day - enrolment_day,
    min_gap   = 365,
    max_inf   = max_inf)

  df <- data.frame(
    subject_id = subject_id,
    time       = sampling_times)

  log_titre_vals <- rep(0, n_samples)

  event_times <- c(0, infection_times)

  n_events <- length(event_times)


  for(j in 1:n_events)
  {
    if(serostatus == 0 & j == 1) next

    start_interval <- event_times[[j]]

    end_interval <- ifelse(j == n_events,
                           sampling_times[[n_samples]],
                           event_times[[j + 1]])

    lower_bound    <- min(sampling_times[sampling_times > start_interval])
    upper_bound    <- max(sampling_times[sampling_times <= end_interval])

    if(j == 1) lower_bound <- sampling_times[[1]]

    interval_times <- sampling_times[sampling_times >= lower_bound &
                                       sampling_times <= upper_bound]

    time_since_event <- interval_times - start_interval

    interval_idx   <- which(sampling_times %in% interval_times)

    sim_titres <- titre_decay_floor(par_alpha = 6,
                                    par_beta  = 2,
                                    par_delta = 0.003,
                                    time      = time_since_event)

    log_titre_vals[interval_idx] <- sim_titres
  }

  df$true_titre <- log_titre_vals

  df
}
