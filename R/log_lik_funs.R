calculate_log_lik <- function(infection_times, sampling_times,
                         perm_rise, temp_rise, temp_decay,
                         measurements, meas_sd)
{
  sim_vals <- simulate_titres_from_enrolment(
    infection_times = infection_times,
    sampling_times  = sampling_times,
    perm_rise       = perm_rise,
    temp_rise       = temp_rise,
    temp_decay      = temp_decay, baseline = measurements[[1]])

  ll_val <- dnorm(x = measurements, mean = sim_vals, sd = meas_sd,
                  log = TRUE) |>
    sum()
}


