#' Simulate Antibody Titre Trajectories for an Individual
#'
#' Generates a simulated time course of true and measured antibody titres
#' for an individual, given sampling times, infection/vaccination events,
#' and baseline serostatus. The function applies a rise-and-decay model
#' with measurement error and limit of detection.
#'
#' @param sampling_times Numeric vector of times at which titres are sampled.
#' @param serostatus Integer (0 or 1) indicating baseline serostatus at
#'  enrolment (0 = seronegative, 1 = seropositive).
#' @param treatment_group Integer (0 or 1) indicating treatment group.
#'  (0 = Placebo, 1 = Vaccine).
#' @param subject_id Identifier for the subject (numeric or character).
#' @param meas_sd A non-negative real value that corresponds to the standard
#'   deviation of a normal measurement model of log titres.
#' @param infection_times Numeric vector of times of infection events relative
#'   to enrolment. Defaults to `NULL` if no infections.
#' @param vac_times Numeric vector of times of vaccination events relative to
#'   enrolment. Defaults to `NULL` if no vaccinations.

#'
#' @returns A `data.frame` with columns:
#' \describe{
#'   \item{subject_id}{Subject identifier.}
#'   \item{time}{Sampling time.}
#'   \item{true}{Simulated true log-titre at each time.}
#'   \item{meas}{Observed log-titre with measurement error and LOD applied.}
#' }
#' @export
#'
#' @examples
#'   subject_id      <- 1
#'   sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#'   serostatus      <- 0 # Seronegative
#'   treatment_group <- 0 # Placebo
#'   simulate_titre_trajectory(sampling_times,
#'                             serostatus,
#'                             treatment_group,
#'                             subject_id      = subject_id,
#'                             infection_times = 1047.108)
simulate_titre_trajectory <- function(sampling_times, serostatus,
                                      treatment_group,
                                      subject_id,
                                      meas_sd,
                                      infection_times = NULL,
                                      vac_times       = NULL)
{
  n_samples     <- length(sampling_times)

  df <- data.frame(
    subject_id = subject_id,
    time       = sampling_times)

  log_titre_vals <- rep(0, n_samples)

  event_times <- c(0, infection_times)

  n_events <- length(event_times)

  # First event corresponds to enrolment
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

    times_within_interval <- sampling_times[sampling_times >= lower_bound &
                                                 sampling_times <= upper_bound]

    # We add `end_interval` to the estimate the last titre right before the next
    #  event
    interval_times <- c(times_within_interval, end_interval)

    time_since_event <- interval_times - start_interval

    interval_idx   <- which(sampling_times %in% times_within_interval)

    perm_rise <- 6 # Permanent rise

    temp_rise_seroneg <- 2 # Temporary rise for seronegative individuals

    decay_rate <- 0.003

    delta_seropos <- ifelse(serostatus == 1, last_titre_prev_inf - perm_rise, 0)

    temp_rise <- temp_rise_seroneg + delta_seropos

    sim_titres <- titre_decay_floor(par_alpha = perm_rise,
                                    par_beta  = temp_rise,
                                    par_delta = decay_rate,
                                    time      = time_since_event)

    log_titre_vals[interval_idx] <- sim_titres[-length(sim_titres)]

    # Last titre of the previous infection
    last_titre_prev_inf <- sim_titres[length(sim_titres)]

    # Any event (vac or inf) will render an individual seropositive
    serostatus <- 1
  }

  df$true <- log_titre_vals

  df$meas <- measurement_model(log_titre_vals,
                               measurement_error = meas_sd,
                               LOD = 1)

  df
}

# Gaussian measurement errors
measurement_model <- function(true_titre, measurement_error, LOD)
{
  obs_titre <- true_titre + stats::rnorm(n    = length(true_titre),
                                         mean = 0,
                                         sd   = measurement_error)

  obs_titre <- ifelse(obs_titre < LOD, 0, obs_titre)

  obs_titre
}
