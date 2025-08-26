#' Simulate Antibody Titre Trajectories for an Individual
#'
#' Generates a simulated time course of true and measured antibody titres
#' for an individual, given sampling times, infection/vaccination events,
#' and baseline serostatus. The function applies a rise-and-decay model
#' with measurement error and limit of detection.
#'
#' @param sampling_times Numeric vector of sampling times (in days) relative to
#'  the start of the study.
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
#'   sampling_times  <- c(182, 210, 266, 294, 434, 643, 980, 1347, 1740)
#'simulate_titre_trajectory(sampling_times   = sampling_times,
#'                          treatment_group  = 0, # Placebo
#'                          subject_id       = 1,
#'                          infection_times  = c(300, 900),
#'                          meas_sd          = 0.3)
simulate_titre_trajectory <- function(sampling_times,
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

  # 0 corresponds to the date we start tracking the individual.
  # This can be the enrolment date in a study or the date she became susceptible.
  event_times <- c(0, infection_times)

  n_events <- length(event_times)

  serostatus <- 0

  if(length(event_times) > 1)
  {
    for(j in 2:n_events)
    {
      start_interval <- event_times[[j]]

      end_interval <- ifelse(j == n_events,
                             sampling_times[[n_samples]] + 1,
                             event_times[[j + 1]])

      times_within_interval <- sampling_times[sampling_times >= start_interval &
                                                sampling_times <= end_interval]

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


#' Simulate antibody titre trajectories from birth
#'
#' @inheritParams simulate_titre_trajectory
#' @param inf_times_sbs Numeric vector of infection times (in days) relative to
#'  day the individual became susceptible.
#' @param age Numeric. Age of the individual at the time of the first sample
#'  (in years).
#'
#' @returns A data frame containing simulated titre values across the specified
#'  sampling times.
#' @export
#'
#' @examples
#' simulate_titres_seropositive(
#'   inf_times_sbs   = c(1059.646, 2743.543, 4152.091),
#'   sampling_times  = c(41, 74, 122, 157, 290, 468, 787, 1192, 1553),
#'   age             = 10,
#'   subject_id      = 3,
#'   treatment_group = 0,
#'   meas_sd         = 0)
simulate_titres_seropositive <- function(inf_times_sbs, sampling_times, age,
                                        subject_id, treatment_group,
                                        meas_sd) {
  # We assume that the first blood drawn was taken on the individual's birthday
  sampling_times_rel_to_last_birthday <- sampling_times - min(sampling_times)

  # Duration of protection of maternal antibodies
  d_mAB <- 1

  # Sampling times relative to the day the individual became susceptible
  sampling_times_sbs <- sampling_times_rel_to_last_birthday  +
    365 * (age - d_mAB)

  sim_titre_df <- simulate_titre_trajectory(
    sampling_times  = sampling_times_sbs,
    treatment_group = treatment_group,
    subject_id      = subject_id,
    infection_times = inf_times_sbs,
    meas_sd         = meas_sd)

  sim_titre_df$time <- sampling_times

  sim_titre_df
}
